#partition=p_medium
#nprocshared=8
#mem=2300MB

import argparse
import os
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from scripts_for_molecule_database_calculations import update_db, folder_exist, sanitize

import ase.db as db
from ase.calculators.vasp import Vasp
from ase.parallel import parprint, world, barrier


def main(db_id: int, db_dir: str):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        functional = row.get('xc')
        atoms = row.toatoms()
        name = row.get('name')

    parprint(f'outstd of opt calculation for db entry {db_id} with structure: {name} and functional: {functional}')

    barrier()

    calc = Vasp(atoms=atoms,
                txt=f'{name}_{functional}_opt.txt',
                xc=functional,
                system=f'{name}_{functional}_opt',  # Insert the name of the system
                istart=0,  # Wavefunction
                icharg=2,  # Charge: 1-file 2-atom 10-const
                # PREC = Normal           # Specifies the "precision"-mode
                ispin=1,  # Specifies spin polarization (1 no, 2 yes)
                encut=450,  # Energy limit for plane waves in eV (never <300)
                nelmin=6,  # Minimum number of electronic steps
                nelm=500,
                ismear=0,  # Smearing
                sigma=0.001,  # Smearing's Width (10^-3 for molecule)
                ibrion=2,  # Atomic optimization algorithm
                nsw=500,  # Maximum number of ionic steps
                potim=0.25,  # IBRION's scaling factor
                ediff=0.00001,  # Electronic optimization criteria (energy, eV)
                ediffg=-0.01,  # Atomic optimization criteria (negative: forces eV/A, positive: energy eV)
                # npar=4,  # number of bands that are treated in parallel, sqrt of n of cores
                )

    calc.calculate(atoms)

    if world.rank == 0: update_db(db_dir, dict(id=db_id, atoms=calc.atoms, relaxed=True, vibration=False, vib_en=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)
