import argparse
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))
from error_project_san_sebastion import sanitize, folder_exist, update_db
#from error_project_san_sebastion.functional_groups import count_methyls, count_amines, count_iso_carbon, count_neo_carbons

#from ase.io import read
import ase.db as db
from ase import Atoms
from ase.io import write
#from ase.geometry import find_mic
#from ase.visualize import view
#from ase.filters import Filter
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
import numpy as np
#import mofun


def main(db_dir: str, pattern: str, ):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        for row in db_obj.select(selection=pattern):
            atoms: Atoms = row.toatoms()
            name = row.get('name')
            id = row.get('id')
            write(f'POSCAR.{id}_{name}', atoms, format='vasp')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database', help='directory to the database.')
    parser.add_argument('pattern')
#    parser.add_argument('not_pattern')
    args = parser.parse_args()

    main(args.database, args.pattern)