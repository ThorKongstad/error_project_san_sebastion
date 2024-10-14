#partition=p_medium
#nprocshared=12
#mem=4000MB

import argparse
import os
import sys
import pathlib
from subprocess import call

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import update_db, folder_exist, sanitize

import ase.db as db
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import string_keys
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

    string_keys.insert(-1, 'libxc1')
    string_keys.insert(-1, 'libxc2')
    Vasp.xc_defaults['n12'] = dict(gga='LIBXC', libxc1='GGA_X_N12', libxc2='GGA_C_N12')
    Vasp.xc_defaults['mn12l'] = dict(metagga='LIBXC', libxc1='MGGA_X_MN12_L', libxc2='MGGA_C_MN12_L', lasph=True)

    Vasp.xc_defaults['tpss'] = dict(metagga='TPSS', lasph=True)#, algo='A')
    Vasp.xc_defaults['m06l'] = dict(metagga='M06L', lasph=True)#, algo='A')

    if functional in ['tpss', 'm06l']: call(['export', 'VASP_PP_PATH=/home-nas/waaguest/VASP_PP_031024'])

    if all([
       name == 'oxygen'
            ]):
        spin = 2
    else:
        spin = 1

    barrier()

    calc = Vasp(atoms=atoms,
                txt=f'{name}_{functional}_opt.txt',
                xc=functional,
                system=f'{name}_{functional}_opt',  # Insert the name of the system
                istart=0,  # Wavefunction
                icharg=2,  # Charge: 1-file 2-atom 10-const
                # PREC = Normal           # Specifies the "precision"-mode
                ispin=spin,  # Specifies spin polarization (1 no, 2 yes)
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
                lcharg=False,
                lwave=False,
    )

    atoms.calc = calc

    calc.calculate(atoms)

    if world.rank == 0 and calc.converged: update_db(db_dir, dict(id=db_id, atoms=atoms, relaxed=True, vibration=False, vib_en=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('database', help='directory to the database.')
    args = parser.parse_args()

    main(args.data_base_id, args.database)
