import argparse
import os
#import sys
#import pathlib

#sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
#from error_project_san_sebastion import sanitize, folder_exist

#from ase.io import read
import ase.db as db
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
# import numpy as np


def main(backup_database: str, backup_database_id: int, overwrite_database: str, overwrite_database_id: int):
    if not os.path.basename(backup_database) in os.listdir(db_path if len(db_path := os.path.dirname(backup_database)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(backup_database) as db_obj:
        row = db_obj.get(selection=f'id={backup_database_id}')

    # connect to db
    if not os.path.basename(overwrite_database) in os.listdir(db_path if len(db_path := os.path.dirname(overwrite_database)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(overwrite_database) as db_obj:
        # db_id = db_obj.reserve(xc = functional, smiles=smile)
        db_obj.write(id=overwrite_database_id, atoms=row, key_value_pairs=row.key_value_pairs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('backup_database')
    parser.add_argument('backup_database_id', type=int)
    parser.add_argument('overwrite_database')
    parser.add_argument('overwrite_database_id', type=int)
    args = parser.parse_args()

    main(**args.__dict__)
