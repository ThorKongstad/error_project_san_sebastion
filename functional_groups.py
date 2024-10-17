import argparse
import os
import sys
import pathlib
from copy import deepcopy as dcp
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


def sp(x):
    print(x)
    return x


def invert_pos(pos):
    inversion = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])
    inversion_func = lambda pos: inversion.dot(np.array(pos, ndmin=2).transpose()).flatten()
    return tuple(map(inversion_func, pos))


def counter(pattern: mofun.Atoms, structure: mofun.Atoms, not_patterns: mofun.Atoms | tuple[mofun.Atoms, ...] | None = None, first_unique: bool = True, atol: float = 0.05) -> int:
    results: list = mofun.find_pattern_in_structure(structure, pattern, atol=atol)
    if not_patterns is not None:
        if not isinstance(not_patterns, (tuple, list)): not_patterns = (not_patterns,)
        for not_pattern in not_patterns:
            not_results: list = mofun.find_pattern_in_structure(structure, not_pattern, atol=atol)
            not_center_atoms = {res[0] for res in not_results}
            for res in dcp(results):
                if res[0] in not_center_atoms: results.remove(res)
    if first_unique:
        center_atoms = set()
        for res in dcp(results):
            if res[0] in center_atoms: results.remove(res)
            else: center_atoms.add(res[0])
    return len(results)


def count_methyls(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    methyls_count = counter(
        pattern=(methyl_patt := mofun.Atoms(elements="CHHH", positions=(methyl_pos := [[ 0.        ,  0.        ,  0.        ], [ 0.63284602,  0.63284602,  0.63284602], [ 0.63284602, -0.63284602, -0.63284602], [-0.63284602,  0.63284602, -0.63284602]]))),
        structure=mofun_atoms,
        atol=0.2
    )

    alkane_count = counter(
        pattern=mofun.Atoms(elements='CHH', positions=[[-3.18155796e-05,  5.86340743e-01,  6.05774660e-05], [-3.32088576e-05,  1.25011122e+00,  8.80447512e-01], [-6.47463927e-05,  1.25027547e+00, -8.80200664e-01]]),
        structure=mofun_atoms,
        not_patterns=(methyl_patt,),
        atol=0.2
    )

    #patterns = [
    #    mofun.Atoms(elements="CHHH", positions=[[ 0.        ,  0.        ,  0.        ], [ 0.63284602,  0.63284602,  0.63284602], [ 0.63284602, -0.63284602, -0.63284602], [-0.63284602,  0.63284602, -0.63284602]]),
    #    mofun.Atoms(elements="CCCHH", positions=[(0.867, 1.760, 1.589), (2.184, 2.409, 1.163), (3.417, 1.604, 1.580), (2.252, 3.424, 1.591), (2.193, 2.542, 0.067)]),
    #    mofun.Atoms(elements="NCCHH", positions=(pos := [(0.867, 1.760, 1.589), (2.184, 2.409, 1.163), (3.417, 1.604, 1.580), (2.252, 3.424, 1.591), (2.193, 2.542, 0.067)])),
    #    mofun.Atoms(elements="NCCHH", positions=invert_pos(pos))
    #    ]

    return methyls_count + alkane_count


def count_iso_carbon(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)

    alkane_pattern = mofun.Atoms(elements="CHH", positions=[[ 0.        ,  0.        ,  0.        ], [ 0.63284602,  0.63284602,  0.63284602], [ 0.63284602, -0.63284602, -0.63284602]])

    iso_count = counter(
        pattern=mofun.Atoms(elements="CH", positions=[[ 0.        ,  0.        ,  0.        ], [ 0.63284602,  0.63284602,  0.63284602]]),
        structure=mofun_atoms,
        not_patterns=(alkane_pattern,),
        atol=0.2
    )

    return iso_count


def count_neo_carbons(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    neo_carbon = mofun.Atoms(elements='CCCC',positions=(poss := np.array([[17.11237229, 17.11237229, 18.88762771],
           [17.11237229, 18.88762771, 17.11237229],
           [18.88762771, 17.11237229, 17.11237229],
           [18., 18., 18.],
           [18.88762771, 18.88762771, 18.88762771]])))

    neo_carbon_results = mofun.find_pattern_in_structure(mofun_atoms, neo_carbon)
    return len(neo_carbon_results)


def count_amines(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    amine = mofun.Atoms(elements='NHH', positions=(poss := np.array([[ 0.00000000e+00, -9.28668405e-06,  4.03847481e-03],
       [ 0.00000000e+00, -9.45924457e-01, -3.82932255e-01],
       [ 8.19109801e-01,  4.73016872e-01, -3.82953110e-01]])))

    return counter(amine, mofun_atoms)
