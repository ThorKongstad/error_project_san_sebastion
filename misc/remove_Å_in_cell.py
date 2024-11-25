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
#from ase.geometry import find_mic
#from ase.visualize import view
#from ase.filters import Filter
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
import numpy as np
#import mofun


def main(db_id: int, db_dir: str, remove: float = 5):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        functional = row.get('xc')
        atoms: Atoms = row.toatoms()
        name = row.get('name')

    org_height = atoms.cell.cellpar()[2]
    atoms.set_cell([[atoms.cell.cellpar()[0], 0, 0], [0, atoms.cell.cellpar()[1], 0], [0, 0, atoms.cell.cellpar()[2] - 5]], scale_atoms=False)
    for atom in atoms:
        if atom.position[2] >= org_height / 2:
            atom.translate([0, 0, -remove])

    update_db(db_dir, dict(id=db_id, atoms=atoms))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)
