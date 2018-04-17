"""This module is a Python/Cython/C++ implementation of the BUDE force field."""

import pathlib

from .force_field import assign_force_field, BuffForceField
from .calculate_energy import (score_interactions, find_inter_ampal,
                               find_intra_ampal)


FF_PATH = pathlib.Path(__file__).parent / 'force_fields'
FORCE_FIELDS = {}


for _ff in FF_PATH.glob('*.json'):
    _ff_id = _ff.name.split('.')[0]
    FORCE_FIELDS[_ff_id] = BuffForceField(str(_ff), _ff_id)


def get_interaction_energy(ampal_objs, ff=None, assign_ff=True):
    """Calculates the interaction energy between AMPAL objects.

    Parameters
    ----------
    ampal_objs: [AMPAL Object]
        A list of any AMPAL objects with `get_atoms` methods.
    ff: BuffForceField, optional
        The force field to be used for scoring. If no force field is
        provided then the most current version of the BUDE force field
        will be used.
    assign_ff: bool, optional
        If true, then force field assignment on the AMPAL object will be
        will be updated.

    Returns
    -------
    BUFF_score: BUFFScore
        A BUFFScore object with information about each of the interactions and
        the atoms involved.
    """
    if ff is None:
        ff = FORCE_FIELDS['bude_2016v1']
    if assign_ff:
        for ampal_obj in ampal_objs:
            assign_force_field(ampal_obj, ff)
    interactions = find_inter_ampal(ampal_objs, ff.distance_cutoff)
    buff_score = score_interactions(interactions, ff)
    return buff_score


def get_internal_energy(ampal_obj, ff=None, assign_ff=True):
    """Calculates the internal energy of the AMPAL object.

    Parameters
    ----------
    ampal_obj: AMPAL Object
        Any AMPAL object with a `get_atoms` method.
    ff: BuffForceField, optional
        The force field to be used for scoring. If no force field is
        provided then the most current version of the BUDE force field
        will be used.
    assign_ff: bool, optional
        If true, then force field assignment on the AMPAL object will be
        will be updated.

    Returns
    -------
    BUFF_score: BUFFScore
        A BUFFScore object with information about each of the interactions and
        the atoms involved.
    """
    if ff is None:
        ff = FORCE_FIELDS['bude_2016v1']
    if assign_ff:
        assign_force_field(ampal_obj, ff)
    interactions = find_intra_ampal(ampal_obj, ff.distance_cutoff)
    buff_score = score_interactions(interactions, ff)
    return buff_score


del pathlib
