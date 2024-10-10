import argparse
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import sanitize, folder_exist

from ase.io import read
import ase.db as db
from ase import Atoms
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
# import numpy as np
import mofun


def main(db_id: int, db_dir: str):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        functional = row.get('xc')
        atoms: Atoms = row.toatoms()
        name = row.get('name')

    atoms.set_positions(atoms.get_positions(wrap=True, pretty_translation=True))

    mofun_atoms = mofun.Atoms.from_ase_atoms(atoms)

    methyl = mofun.Atoms(elements="CHHH", positions=[(0.867, 1.760, 1.589), (0, 2.36, 1.276), (0.816, 1.646, 2.683), (0.757, 0.758, 1.148)])

    methyl_results = mofun.find_pattern_in_structure(mofun_atoms, methyl)
    print(methyl_results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)
