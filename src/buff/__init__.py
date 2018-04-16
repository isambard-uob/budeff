"""This module is a Python/Cython/C++ implementation of the BUDE force field."""

import pathlib

from .force_field import BuffForceField
from .calculate_energy import (PyAtomData, score_interactions,
                               get_within_ff_cutoff, find_inter_ampal,
                               find_intra_ampal)


FF_PATH = pathlib.Path(__file__).parent / 'force_fields'
FORCE_FIELDS = {}

for _ff in FF_PATH.glob('*.json'):
    _ff_id = _ff.name.split('.')[0]
    FORCE_FIELDS[_ff_id] = BuffForceField(str(_ff), _ff_id)

del pathlib
