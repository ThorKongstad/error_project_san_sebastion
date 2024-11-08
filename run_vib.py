#partition=p_medium
#nprocshared=12
#mem=4000MB

import argparse
import os
import time
import sys
import pathlib
from subprocess import call
from typing import NoReturn

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from scripts_for_molecule_database_calculations import update_db, folder_exist, sanitize, ends_with

import ase.db as db
from ase.vibrations import Vibrations
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import string_keys
from ase.parallel import parprint, world, barrier



def file_dont_exist(file_name: str, path: str = '.', rm_flags='', return_path: bool = False) -> NoReturn | str:
    if file_name in os.listdir(path):
        call(['rm', f'{rm_flags}', f"'{ends_with(path, '/')}{file_name}'"])
    if return_path: return f'{ends_with(path, "/")}{file_name}'


def clean_old_files(functional_folder, file_name):
    os.getcwd()
    if file_name in os.listdir(functional_folder):
        folder_exist(f'old_vibs', path=functional_folder)
        call(['mv', '-f', f"{ends_with(functional_folder, '/')}{file_name}", f"{file_dont_exist(file_name, f'{functional_folder}/old_vibs', return_path=True)}"])
    if (old_folder := file_name.replace('.txt', '')) in os.listdir(functional_folder):
        folder_exist(f'old_vibs', path=functional_folder)
        call(['mv', '-f', f"{ends_with(functional_folder, '/')}{old_folder}", f"{file_dont_exist(old_folder, f'{functional_folder}/old_vibs', rm_flags='-r', return_path=True)}"])


def main(db_id: int, clean_old: bool = True, db_dir: str = 'molreact.db'):
    # read from  database
    if not os.path.basename(db_dir) in os.listdir(db_path if len(db_path := os.path.dirname(db_dir)) > 0 else '.'): raise FileNotFoundError("Can't find database")
    with db.connect(db_dir) as db_obj:
        row = db_obj.get(selection=f'id={db_id}')
        if row.get('relaxed'): atoms = row.toatoms()
        else: raise f"atoms at row id: {db_id} haven't been relaxed."
        name = row.get('name')
        functional = row.get('xc')
        data: dict = row.get('data')

    parprint(f'outstd of vib calculation for db entry {db_id} with structure: {name} and functional: {functional}')

    string_keys.insert(-1, 'libxc1')
    string_keys.insert(-1, 'libxc2')
    Vasp.xc_defaults['n12'] = dict(gga='LIBXC', libxc1='GGA_X_N12', libxc2='GGA_C_N12')
    Vasp.xc_defaults['mn12l'] = dict(metagga='LIBXC', libxc1='MGGA_X_MN12_L', libxc2='MGGA_C_MN12_L', lasph=True)

    Vasp.xc_defaults['tpss'] = dict(metagga='TPSS', lasph=True)  # , algo='A')
    Vasp.xc_defaults['m06l'] = dict(metagga='M06L', lasph=True)  # , algo='A')

    if functional in ['tpss', 'm06l']: os.environ[
        'VASP_PP_PATH'] = '/home-nas/waaguest/VASP_PP_031024'  # call(['export', 'VASP_PP_PATH=/home-nas/waaguest/VASP_PP_031024'])

    if all([
        name == 'oxygen',
        name == 'nitric-oxide',
            ]):
        spin = 2
    else:
        spin = 1

    barrier()

    calc = Vasp(atoms=atoms,
                txt=f'{name}_id{db_id}_{functional}_vib.txt',
                xc=functional,
                system=f'{name}_id{db_id}_{functional}_vib',  # Insert the name of the system

                prec='Normal',
                icharg=2,
                ispin=spin,  # Specifies spin polarization (1 no, 2 yes)
                encut=450,  # Energy limit for plane waves in eV (never <300)
                nelm=500,
                nelmin=6,  # Minimum number of electronic steps
                ismear=0,  # Smearing
                sigma=0.001,  # Smearing's Width (10^-3 for molecule)
                ibrion=5,
                nsw=500,  # Maximum number of ionic steps
                potim=0.02,  # IBRION's scaling factor
                ediff=0.00001,  # Electronic optimization criteria (energy, eV)
                ediffg=-0.01,  # Atomic optimization criteria (negative: forces eV/A, positive: energy eV)
                lwave=False,
                lcharg=False,
                nfree=2,
                isym=-1,
                algo='A',
                istart=0,  # Wavefunction
                )

    atoms.calc = calc

    calc.calculate(atoms)
    freq, ifreq = calc.read_vib_freq()
    zpe = sum(freq)


    if world.rank == 0:
        data.update(dict(frequencies=freq, i_frequencies=ifreq))
    #    vib.summary(log=f'{functional_folder}/{file_name.replace("vib", "vib_en")}')

    #    with open(f'{functional_folder}/{file_name.replace("vib", "vib_en")}', 'r') as fil:
    #        energy_string = fil.read()

        # saving vib data
        update_db(db_dir, db_update_args=dict(
            id=db_id,
            vibration=True,
            zpe=zpe,
            data=data
        ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_base_id', type=int)
    parser.add_argument('-db', '--database', help='directory to the database, if not stated will look for molreact.db in pwd.', default='molreact.db')
    args = parser.parse_args()

    main(args.data_base_id, db_dir=args.database)
