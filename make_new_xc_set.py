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


def main(model_xc: str, new_xc: str, db_dir: str):
    with db.connect(db_dir) as db_obj:
        row_iter = db_obj.select(selection=f'xc={model_xc}',)
        for row in row_iter:
            atoms = row.toatoms().copy()
            name = row.get('name')
            db_obj.write(atoms=atoms, xc=new_xc, name=name, relaxed=False, vibration=False, )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('db_dir',)
    parser.add_argument('model_xc')
    parser.add_argument('new_xc')
    args = parser.parse_args()

    main(**args.__dict__)
