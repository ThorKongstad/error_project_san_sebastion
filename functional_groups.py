from copy import deepcopy as dcp
import numpy as np
from ase import Atoms
import mofun


def sp(x):
    print(x)
    return x


def invert_pos(pos):
    inversion = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])
    inversion_func = lambda pos: inversion.dot(np.array(pos, ndmin=2).transpose()).flatten()
    return tuple(map(inversion_func, pos))


def big_pattern_catcher(func):
    def wrap(*args, **kwargs):
        try: return func(*args, **kwargs)
        except Exception: return 0
    return wrap


@ big_pattern_catcher
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
        atol=0.5
    )

    alkane_count = counter(
        pattern=mofun.Atoms(elements='CHH', positions=[[-3.18155796e-05,  5.86340743e-01,  6.05774660e-05], [-3.32088576e-05,  1.25011122e+00,  8.80447512e-01], [-6.47463927e-05,  1.25027547e+00, -8.80200664e-01]]),
        structure=mofun_atoms,
        not_patterns=(methyl_patt,),
        atol=0.5
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
        atol=0.5
    )

    return iso_count


def count_neo_carbons(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    neo_carbon = mofun.Atoms(elements='CCCCC', positions=(poss := np.array([[0.0, 0.0, 0.0, ],
                                                        [ 0.8876277091687248, 0.8876277091687248, 0.8876277091687248, ],
                                                        [ -0.887627709168727, -0.887627709168727, 0.8876277091687248, ],
                                                        [ -0.887627709168727, 0.8876277091687248, -0.887627709168727, ],
                                                        [ 0.8876277091687248, -0.887627709168727, -0.887627709168727]])))

    neo_carbon_results = mofun.find_pattern_in_structure(mofun_atoms, neo_carbon)
    return len(neo_carbon_results)


def count_amines(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    amine = mofun.Atoms(elements='NHH', positions=(poss := np.array([[ 0.00000000e+00, -9.28668405e-06,  4.03847481e-03],
       [ 0.00000000e+00, -9.45924457e-01, -3.82932255e-01],
       [ 8.19109801e-01,  4.73016872e-01, -3.82953110e-01]])))

    return counter(amine, mofun_atoms)


def count_nitro(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    nitrate = mofun.Atoms(elements='NOOO',
                          positions=(poss := np.array([[0.011908105419858003, 0.6033905189532001, 0.0, ],
                                                       [-1.1972268438180045, 0.5643155515465921, 0.0, ],
                                                       [0.746479277814504, 1.5547357964356223, 0.0, ],
                                                       [0.6911508880881579, -0.6097007037117663, 0.0]])))
    nitro = mofun.Atoms(elements='NOO', positions=(poss := np.array([[ -0.09, -0.006, 0],
       [-0.633, 1.086, 0],
       [-0.6622409883469782, -1.084061319218641, 0.0]])))

    return counter(nitro, mofun_atoms, not_patterns=nitrate, atol=0.5)


def count_nitrate(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    nitrate = mofun.Atoms(elements='NOOO', positions=(poss := np.array([[0.011908105419858003, 0.6033905189532001, 0.0, ],
                                                                        [-1.1972268438180045, 0.5643155515465921, 0.0, ],
                                                                        [0.746479277814504, 1.5547357964356223, 0.0, ],
                                                                        [0.6911508880881579, -0.6097007037117663, 0.0]])))

    return counter(nitrate, mofun_atoms, atol=0.2)


def count_hydroxylamine(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    hydroxylamine = mofun.Atoms(elements='NHO', positions=(poss := np.array([[-0.04137370584640099, 0.017367187728555, 0.20813707658508657,],
                                                                        [-0.07674818036593378, -1.8822985589231374, 0.3622425058205117,],
                                                                        [-0.067365431436752, -1.2842511498652631, -0.39161222155933206]])))

    return counter(hydroxylamine, mofun_atoms, atol=0.2)


def count_hydrazine(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    hydrazine = mofun.Atoms(elements='NN', positions=(poss := np.array([[-0.05104239599528282, 0.7165023983018459, -0.0017242674984179196],
                                                                        [0.051042395995282805, -0.7165023983018488, -0.0017242674984179196]])))

    return counter(hydrazine, mofun_atoms, atol=0.2)


def count_nitrile(atoms: Atoms | mofun.Atoms) -> int:
    mofun_atoms = atoms if isinstance(atoms, mofun.Atoms) else mofun.Atoms.from_ase_atoms(atoms)
    nitrile = mofun.Atoms(elements='NC', positions=(poss := np.array([[0.0003618272271315, 0.00064023551643, -1.4296814167539034],
                                                                      [9.49608749745e-05, 0.0001270950249525, -0.27591085098963575]])))
    Nitrogen_hydrogen = mofun.Atoms(elements='NH', positions=(poss := np.array([[0.00000000e+00, -9.28668405e-06,  4.03847481e-03],
       [0.00000000e+00, -9.45924457e-01, -3.82932255e-01]])))

    return counter(nitrile, mofun_atoms, not_patterns=Nitrogen_hydrogen, atol=0.2)

