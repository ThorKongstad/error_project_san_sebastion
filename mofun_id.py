import argparse
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import sanitize, folder_exist
from error_project_san_sebastion.functional_groups import count_methyls, count_amines, count_iso_carbon, count_neo_carbons, count_nitro, count_nitrate, count_hydroxylamine, count_hydrazine, count_nitrile

from ase.io import read
import ase.db as db
from ase import Atoms
from ase.geometry import find_mic
#from ase.filters import Filter
#from ase.data.pubchem import pubchem_atoms_search
#from ase.build import molecule
import numpy as np
import mofun


def get_repeated_representation(atoms: Atoms, selection_dis: float = 7) -> Atoms:
    # a try of repeating the structure and then only taking the atom that apears in the middle of the repeat
    atoms_repeated = atoms.repeat((2, 2, 2))
    cell_corner = sum(atoms.cell.array)
    sphere_filter = lambda atom: np.linalg.norm(np.array(atom.position) - cell_corner) < selection_dis
    #rep_atoms = Filter(atoms_repeated, mask=list(map(sphere_filter, atoms_repeated)))
    return atoms_repeated[list(map(sphere_filter, atoms_repeated))]


def main(db_id: int, db_dir: str):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        functional = row.get('xc')
        atoms: Atoms = row.toatoms()
        name = row.get('name')

#    atoms.set_positions(atoms.get_positions(wrap=True, pretty_translation=True))

#    repr_atoms = get_repeated_representation(atoms)
#    assert len(repr_atoms) == len(atoms)

    atoms.set_positions(find_mic(atoms.get_positions(), atoms.cell)[0])

    mofun_atoms = mofun.Atoms.from_ase_atoms(atoms)

    print(f'{name}: There are {count_methyls(mofun_atoms)} methyl groups, {count_iso_carbon(mofun_atoms)} iso-carbons, {count_neo_carbons(mofun_atoms)} neo-carbons.')

    print(f' there are {count_amines(mofun_atoms)} amines, {count_nitro(mofun_atoms)} nitro groups, {count_nitrate(mofun_atoms)} nitrates, {count_hydroxylamine(mofun_atoms)}, hydroxylamines, {count_hydrazine(mofun_atoms)} hydrazines, and {count_nitrile(mofun_atoms)} nitrile')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)
