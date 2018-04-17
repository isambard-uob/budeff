"""Contains code for BUFF force field objects."""

import json
import pathlib
import warnings

from .calculate_energy import PyAtomData


def assign_force_field(ampal_obj, ff):
    """Assigns force field parameters to Atoms in the AMPAL object.

    Parameters
    ----------
    ampal_obj : AMPAL Object
        Any AMPAL object with a `get_atoms` method.
    ff: BuffForceField
        The force field to be used for scoring.
    """
    if hasattr(ampal_obj, 'ligands'):
        atoms = ampal_obj.get_atoms(ligands=True, inc_alt_states=True)
    else:
        atoms = ampal_obj.get_atoms(inc_alt_states=True)
    for atom in atoms:
        w_str = None
        a_ff_id = None
        if atom.element == 'H':
            continue
        elif atom.parent.mol_code.upper() in ff:
            if atom.res_label.upper() in ff[atom.parent.mol_code]:
                a_ff_id = (atom.parent.mol_code.upper(),
                           atom.res_label.upper())
            elif atom.res_label.upper() in ff['WLD']:
                a_ff_id = ('WLD', atom.res_label.upper())
            else:
                w_str = ('{} atom is not parameterised in the selected '
                         'force field for {} residues, this will be '
                         'ignored.').format(
                             atom.res_label, atom.parent.mol_code)
        elif atom.res_label.upper() in ff['WLD']:
            a_ff_id = ('WLD', atom.res_label.upper())
        else:
            w_str = ('{} ({}) atom is not parameterised in the selected'
                     ' residue force field.').format(
                         atom.res_label, atom.parent.mol_code)
        if w_str:
            warnings.warn(w_str, NotParameterisedWarning)
        atom.tags['_buff_ff_id'] = a_ff_id
    return


class NotParameterisedWarning(RuntimeWarning):
    pass


class ForceFieldParameterError(Exception):
    pass


class BuffForceField(dict):
    """A wrapper around a BUFF force field.

    Parameters
    ----------
    force_field : string
        Name of force field to be loaded.
    auto_update_params : bool, optional

    Attributes
    ----------
    """

    def __init__(self, force_field_path, force_field_name, auto_update_params=False):
        with open(force_field_path, 'r') as inf:
            in_d = json.loads(inf.read())
        super().__init__(in_d)
        self._parameter_struct_dict = None
        self._old_hash = None
        self._defined_dist_cutoff = None
        self.path = force_field_path
        self.name = force_field_name
        self.auto_update_f_params = auto_update_params

    def __repr__(self):
        return "<BUFF Force Field Object: {}>".format(self.name)

    @property
    def max_radius_and_npnp(self):
        """Maximum radius and non-polar non-polar distance in the force field."""
        return self.find_max_rad_npnp()

    @property
    def distance_cutoff(self):
        """Distance cut off for interactions within the force field."""
        if self._defined_dist_cutoff is None:
            return self._calc_distance_cutoff()
        else:
            return self._defined_dist_cutoff

    @distance_cutoff.setter
    def distance_cutoff(self, cutoff):
        self._defined_dist_cutoff = cutoff
        return

    def _calc_distance_cutoff(self):
        rad, npnp = self.find_max_rad_npnp()
        return (rad * 2) + npnp

    def find_max_rad_npnp(self):
        """Finds the maximum radius and npnp in the force field.

        Returns
        -------
        (max_rad, max_npnp): (float, float)
            Maximum radius and npnp distance in the loaded force field.
        """
        max_rad = 0
        max_npnp = 0
        for res, _ in self.items():
            if res != 'KEY':
                for _, ff_params in self[res].items():
                    if max_rad < ff_params[1]:
                        max_rad = ff_params[1]
                    if max_npnp < ff_params[4]:
                        max_npnp = ff_params[4]
        return max_rad, max_npnp

    @property
    def parameter_struct_dict(self):
        """Dictionary containing PyAtomData structs for the force field."""
        if self._parameter_struct_dict is None:
            self._parameter_struct_dict = self._make_ff_params_dict()
        elif self.auto_update_f_params:
            new_hash = hash(
                tuple([tuple(item)
                       for sublist in self.values()
                       for item in sublist.values()]))
            if self._old_hash != new_hash:
                self._parameter_struct_dict = self._make_ff_params_dict()
                self._old_hash = new_hash
        return self._parameter_struct_dict

    def _make_ff_params_dict(self):
        """Makes a dictionary containing PyAtomData for the force field.

        Returns
        -------
        ff_params_struct_dict: dict
            Dictionary containing PyAtomData structs for the force field
            parameters for each atom in the force field.
        """
        try:
            ff_params_struct_dict = {}
            for res in self.keys():
                if res == 'KEY':
                    continue
                if res not in ff_params_struct_dict:
                    ff_params_struct_dict[res] = {}
                for atom, params in self[res].items():
                    ff_params_struct_dict[res][atom] = PyAtomData(
                        atom.encode(), params[0].encode(), *params[1:])
        except TypeError:
            raise ForceFieldParameterError(
                'Badly formatted force field parameters: {}'.format(params))
        return ff_params_struct_dict


__author__ = "Christopher W. Wood"
