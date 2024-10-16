import argparse
import os
import sys
import pathlib

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import sanitize, folder_exist

from ase.io import read
import ase.db as db
from ase import Atoms
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
# import numpy as np
import mofun


def invert_pos(pos):
    inversion = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])

    inversion_func = lambda pos: inversion * np.array(pos)
    return tuple(map(inversion_func, pos))


def counter(pattern, structure) -> int:
    return len(mofun.find_pattern_in_structure(structure, pattern))


def count_methyls(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    patterns = [
        mofun.Atoms(elements="CHHH", positions=[(0.867, 1.760, 1.589), (0, 2.36, 1.276), (0.816, 1.646, 2.683),(0.757, 0.758, 1.148)]),
        mofun.Atoms(elements="CCCHH", positions=[(0.867, 1.760, 1.589), (2.184, 2.409, 1.163), (3.417, 1.604, 1.580), (2.252, 3.424, 1.591), (2.193, 2.542, 0.067)]),
        mofun.Atoms(elements="NCCHH", positions=(pos := [(0.867, 1.760, 1.589), (2.184, 2.409, 1.163), (3.417, 1.604, 1.580), (2.252, 3.424, 1.591), (2.193, 2.542, 0.067)])),
        mofun.Atoms(elements="NCCHH", positions=invert_pos(pos))
        ]

    return sum(counter(patt, mofun_atoms) for patt in patterns)


def count_iso_carbon(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    iso_carbon = mofun.Atoms('CCCH', positions=(pos := [[7.32079917, 1.97808864, 1.09854328],
       [7.64738407, 3.3509043 , 1.69764671],
       [5.96657801, 1.4319005 , 1.5814139 ],
       [8.44228529, 0.97911226, 1.40304422],
       [7.25888048, 2.09304343, 0.        ]]))

    iso_carbon_inverted = mofun.Atoms('CCCH', positions=invert_pos(pos))

    iso_carbon_inverted_results = mofun.find_pattern_in_structure(mofun_atoms, iso_carbon_inverted)
    iso_carbon_results = mofun.find_pattern_in_structure(mofun_atoms, iso_carbon)
    return len(iso_carbon_results) + len(iso_carbon_inverted_results)


def count_neo_carbons(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    neo_carbon = mofun.Atoms('CCCC',positions=(poss := np.array([[17.11237229, 17.11237229, 18.88762771],
           [17.11237229, 18.88762771, 17.11237229],
           [18.88762771, 17.11237229, 17.11237229],
           [18., 18., 18.],
           [18.88762771, 18.88762771, 18.88762771]])))

    neo_carbon_results = mofun.find_pattern_in_structure(mofun_atoms, neo_carbon)
    return len(neo_carbon_results)


def count_amines(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    amine = mofun.Atoms('CNHH', positions=(poss := np.array([[1.33327003, 19.571581  , 0.],
                                                             [ 2.552727  ,  0.38599328,  0.        ],
                                                             [ 2.57050018,  0.99844134,  0.8183404 ],
                                                             [ 2.57050018,  0.99844134, 19.1816596 ]])))

    amine_results = mofun.find_pattern_in_structure(mofun_atoms, amine)
    return len(amine_results)
