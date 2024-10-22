import argparse
import os
import sys
import pathlib

#sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
#from error_project_san_sebastion import sanitize, folder_exist
#from error_project_san_sebastion.functional_groups import count_methyls, count_amines, count_iso_carbon, count_neo_carbons

#from ase.io import read
import ase.db as db
from ase import Atoms
from ase.geometry import find_mic
from ase.visualize import view
#from ase.filters import Filter
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
import numpy as np
import mofun


def main(db_id: int, db_dir: str):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        functional = row.get('xc')
        atoms: Atoms = row.toatoms()
        name = row.get('name')

    atoms.set_positions(find_mic(atoms.get_positions(), atoms.cell)[0])

    view(atoms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)