import argparse
#import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import sanitize, folder_exist

from ase.io import read
import ase.db as db
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
# import numpy as np


def main(struc_path: str, name: str, functional: str, db_dir: str, overwrite=None):
    # create ase mol
    #if 'traj' == struc_path[-4:]: atoms = read(struc_path)
    #else: atoms = read(struc_path, format='vasp')

    atoms = read(struc_path)

    atoms.pbc = True

    # connect to db

    with db.connect(db_dir) as db_obj:
        # db_id = db_obj.reserve(xc = functional, smiles=smile)
        db_obj.write(id=overwrite, atoms=atoms, xc=functional, name=name, relaxed=False, vibration=False,)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('structure')
    parser.add_argument('smiles_str')
    parser.add_argument('functional', help='str denoting what fucntional to calculate with')
    parser.add_argument('database', help='name or directory for the database.',)
    parser.add_argument('-overwrite', type=int, default=None)
    args = parser.parse_args()

    main(struc_path=args.structure, name=args.smiles_str, functional=args.functional, db_dir=args.database)
