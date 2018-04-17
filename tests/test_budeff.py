"""Contains tests for the BUFF module."""

import pathlib
import unittest
import warnings

import ampal
import budeff
import numpy

warnings.filterwarnings("ignore")

TEST_FILE_FOLDER = pathlib.Path(__file__).parent / 'testing_files'


class ForceFieldTestCase(unittest.TestCase):
    """Tests for budeff.BuffForceField."""

    def setUp(self):
        self.ff = budeff.FORCE_FIELDS['bude_2016v1']
        pdb_path = TEST_FILE_FOLDER / '3qy1.pdb'
        self.pdb = ampal.load_pdb(str(pdb_path))

    def test_default_ff(self):
        """Tests if all the parameters are correctly loaded."""
        for res, atoms in self.ff.items():
            if res != "KEY":
                for atom, params in atoms.items():
                    self.assertEqual(len(params), 8)
                    self.assertTrue(isinstance(atom, str))
                    self.assertTrue(isinstance(params[0], str))
                    self.assertTrue(
                        all([isinstance(x, (float, int)) for x in params[1:]])
                    )

    def test_parameterisation_pdb(self):
        """Checks that all atoms in PDB are parameterised."""
        budeff.assign_force_field(self.pdb, self.ff)
        for atom in self.pdb.get_atoms(inc_alt_states=True):
            if atom.element != 'H':
                if atom.parent.mol_code != 'HOH':
                    self.assertTrue(atom.tags['_buff_ff_id'] is not None)


class InteractionsTestCase(unittest.TestCase):
    """Tests for budeff.find_buff_interactions."""

    def setUp(self):
        self.ff = budeff.FORCE_FIELDS['bude_2015v1']
        # Uses old force field for continuity
        pdb_path = TEST_FILE_FOLDER / '3qy1.pdb'
        self.pdb = ampal.load_pdb(str(pdb_path))
        budeff.assign_force_field(self.pdb, self.ff)

    def test_basic_intra_format(self):
        """Tests that the interaction tuples are correctly formatted."""
        buff_interactions = budeff.find_intra_ampal(
            self.pdb[0], self.ff.distance_cutoff)
        for a, b in buff_interactions:
            self.assertTrue(isinstance(a, ampal.Atom))
            self.assertTrue(isinstance(b, ampal.Atom))
            self.assertTrue(a != b)

    def test_basic_inter_format(self):
        """Tests that the interaction tuples are correctly formatted."""
        buff_interactions = budeff.find_inter_ampal(
            self.pdb, self.ff.distance_cutoff)
        for a, b in buff_interactions:
            self.assertTrue(isinstance(a, ampal.Atom))
            self.assertTrue(isinstance(b, ampal.Atom))
            self.assertTrue(a != b)

    def test_interaction_energy(self):
        """Tests the interaction energy of a reference structure."""
        buff_score = budeff.get_interaction_energy(self.pdb, self.ff)
        self.assertAlmostEqual(buff_score.total_energy, -1005.41, places=2)

    def test_internal_energy(self):
        """Tests the internal energy of a reference structure."""
        buff_score = budeff.get_internal_energy(self.pdb[0], self.ff)
        self.assertAlmostEqual(buff_score.total_energy, -3722.49, places=2)

    def test_inter_score_components(self):
        """Tests that inter_scores is equal to the summed components."""
        buff_score = budeff.get_interaction_energy(self.pdb, self.ff)
        steric = 0
        desolvation = 0
        charge = 0
        for _, score in buff_score.inter_scores:
            steric += score[0]
            desolvation += score[1]
            charge += score[2]
        self.assertEqual(buff_score.steric, steric)
        self.assertEqual(buff_score.desolvation, desolvation)
        self.assertEqual(buff_score.charge, charge)
        self.assertTrue(numpy.isclose(
            buff_score.total_energy, sum([steric, desolvation, charge])))

    def test_inter_score_distances(self):
        """Tests that recorded interactions are within ff cutoff distance."""
        buff_score = budeff.get_interaction_energy(self.pdb, self.ff)
        ff_co = self.ff.distance_cutoff
        for (at_a, at_b), _ in buff_score.inter_scores:
            self.assertTrue(ampal.geometry.distance(at_a, at_b) <= ff_co)


if __name__ == '__main__':
    unittest.main()
