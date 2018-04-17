# Bristol University Docking Engine Force Field (BUDEFF)

BUDEFF is a standalone implementation of the
[BUDE](http://www.bristol.ac.uk/biochemistry/research/bude/) (Bristol University
 Docking Engine) all-atom force field<sup>1,2</sup>. The force field is developed by the [Sessions
group](http://www.bris.ac.uk/biochemistry/people/richard-b-sessions/index.html).

[![CircleCI](https://circleci.com/gh/isambard-uob/budeff/tree/master.svg?style=shield)](https://circleci.com/gh/isambard-uob/budeff/)
[![Python Version](https://img.shields.io/badge/python-3.5%2C%203.6-lightgrey.svg)]()
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/isambard-uob/budeff/blob/master/LICENSE)

## Installation

You can install BUDEFF from pip:

`pip install budeff`

Or from source by downloading/cloning this repository, navigating to the folder
and typing:

`pip install .`

BUDEFF uses Cython, so if you're installing from source make sure you have it
installed.

## Usage

The BUDE force field can be used to calculate energies for any protein structure
that has been loaded into [AMPAL](https://github.com/isambard-uob/ampal/), a
simple framework for representing biomolecular structure. You can load a
structure into AMPAL like so:

```Python
import ampal
structure = ampal.load_pdb('3qy1.pdb')
```

Once the structure is loaded in, you can now run BUFF on the structure. BUFF has
two modes:

1. Internal Energy - In this mode a single AMPAL object is supplied and the
   energy between every pair of atoms is calculated, so long as the atom is parameterised in
   the [force field](
   https://github.com/isambard-uob/budeff/tree/master/src/budeff/force_fields/).
1. Interaction Energy - In this mode a list of AMPAL objects is supplied and the
   energy of atom pairs forming the interaction between these objects is
   measured. For example, if the interaction energy between object A and object
   B, then all atom pairs will contain one atom from A and one from B.

```Python
import budeff
internal_energy = budeff.get_internal_energy(structure)
# OUT: NotParameterisedWarning: O (HOH) atom is not parameterised in the selected residue force field.
# OUT:   warnings.warn(w_str, NotParameterisedWarning)
```

While the score was being calculated, a `NotParameterisedWarning` was raised.
This tells us that the water (HOH) is not parameterised in the force field and
so will be ignored. The BUDE force field has been developed for performing
protein docking, and so only protein and a few common ions are parameterised.

`get_internal_energy` returns a `BuffScore` object:

```Python
print(internal_energy)
# OUT: <BUFF Score -7108.00: 214.37 St | -4343.46 De | -2978.91 Ch>
```

The `BuffScore` contains information on the total energy of the system (-7108.00
in this case) as well as the different components of this score, which are
steric (`214.37 St`), energy of desolvation (`-4343.46 De`) and charged
interactions (`-2978.91 Ch`). Each of these components can be accessed
individually:

```Python

print(internal_energy.total_energy, internal_energy.steric,
      internal_energy.desolvation, internal_energy.charge)
# OUT: -7108.000086377617 214.36602045772776 -4343.460484501997 -2978.905622333365
```

Individual pairwise interactions can be examined. The `inter_scores` attribute
is a list of all the pairwise interactions with non-zero scores that are used to
create the score:

```Python
print(internal_energy.inter_scores[0])
# OUT: ((<Carbon Atom (CA). Coordinates: (15.518, -30.153, -25.207)>,
# OUT:   <Carbon Atom (CB). Coordinates: (17.842, -27.509, -21.862)>),
# OUT:  [0.0, -0.10352520993045879, 0.0])
```

Each element in `inter_scores` contains a pair of atoms which form the
interaction and a list with the different elements of the scoring function in
the order steric, desolvation and charge.

To calculate the interaction energy, use the `get_interaction_energy` function.
This take a list of ampal objects and calculates the interaction energy between
these objects:

```Python
interaction_energy = budeff.get_interaction_energy([structure[0], structure[1]])
print(interaction_energy)
# OUT: NotParameterisedWarning: O (HOH) atom is not parameterised in the selected residue force field.
# OUT:   warnings.warn(w_str, NotParameterisedWarning)
# OUT: <BUFF Score -479.44: 26.19 St | -416.31 De | -89.32 Ch>
```

The score is lower in this case as only the energy between an atom in chain a
and an atom in chain b is considered.

There's lots more functionality in the BUFF module so have a dig around.

## References

1. McIntosh-Smith S. et al. (2012) Benchmarking energy efficiency, power costs
   and carbon emissions on heterogeneous systems. *Comput. J.*, 55, 192–205.
2. McIntosh-Smith S. et al. (2014) High performance in silico virtual drug
   screening on many-core processors. *Int. J. High Perform. Comput. Appl.*, 29,
   119–134.
