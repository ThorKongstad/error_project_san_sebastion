import argparse
import os.path
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import update_db
from ase.io import read
import ase.db as db
from ase import Atoms


def main(db_index: int, txt_dir: str, db_dir: str, relaxed: bool = False, new_constraint: bool = False):
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_index}')
        functional = row.get('xc')
        old_atoms = row.toatoms()
        name = row.get('name')
    if not new_constraint:
        constraints = old_atoms.constraints
    if os.path.basename(txt_dir) == 'OUTCAR': updated_atoms: Atoms = read(txt_dir, index=-1, format='vasp-out')
    elif os.path.basename(txt_dir) == 'CONTCAR': updated_atoms: Atoms = read(txt_dir, format='vasp')
    else: updated_atoms: Atoms = read(txt_dir, index=-1)
    if not new_constraint: updated_atoms.set_constraint(constraints)
    assert old_atoms.get_chemical_formula() == updated_atoms.get_chemical_formula()
    update_db(db_dir,
              dict(
                  id=db_index,
                  atoms=updated_atoms,
                  relaxed=relaxed,
                  vibration=False,
                  vib_en=False
              ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('db_index', type=int)
    parser.add_argument('txt_directory', help='Path to the gpaw output.')
    parser.add_argument('db_directory', help='Path to the database')
    parser.add_argument('-relaxed', '--relaxed', action='store_true', default=False)
    parser.add_argument('-newCon', '--new_constraint', action='store_true', default=False)
    args = parser.parse_args()

    main(args.db_index, args.txt_directory, args.db_directory, args.relaxed, args.new_constraint)
