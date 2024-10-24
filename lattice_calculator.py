#partition=p_medium
#nprocshared=24
#mem=4000MB

import argparse
import os
from typing import Tuple, Sequence, NoReturn, Optional, Callable, Never
#import numpy as np
from ase.data import reference_states,chemical_symbols
from ase.build import bulk, graphene#, fcc100, bcc100,hcp0001,
import ase.db as db
from ase.parallel import parprint, world
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import string_keys
from ase.calculators.mixing import SumCalculator
#from dftd4.ase import DFTD4
#from gpaw.cluster import Cluster
#from gpaw import GPAW, PW, Davidson
#from gpaw.utilities import h2gpts
#from collections import namedtuple
from scipy.optimize import curve_fit, minimize, OptimizeResult
import csv
#from kplib import get_kpoints
#from pymatgen.io.ase import AseAtomAdaptor
import pathlib
from time import sleep
from random import randint
from dataclasses import field, make_dataclass
#import plotly.graph_objects as go
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_random


string_keys.insert(-1, 'libxc1')
string_keys.insert(-1, 'libxc2')
Vasp.xc_defaults['n12'] = dict(gga='LIBXC', libxc1='GGA_X_N12', libxc2='GGA_C_N12')
Vasp.xc_defaults['mn12l'] = dict(metagga='LIBXC', libxc1='MGGA_X_MN12_L', libxc2='MGGA_C_MN12_L', lasph=True)

Vasp.xc_defaults['tpss'] = dict(metagga='TPSS', lasph=True)#, algo='A')
Vasp.xc_defaults['m06l'] = dict(metagga='M06L', lasph=True)#, algo='A')


@retry(retry=retry_if_exception_type(FileExistsError), stop=stop_after_attempt(5), wait=wait_random(min=1, max=25))
def folder_exist(folder_name: str, path: str = '.') -> None:
    if folder_name not in os.listdir(path): os.mkdir(ends_with(path, '/') + folder_name)


def ends_with(string: str, end_str: str) -> str:
    return string + end_str * (end_str != string[-len(end_str):0])


def sanitize(unclean_str: str) -> str:
    for ch in ['!', '*', '?', '{', '[', '(', ')', ']', '}', "'", '"']: unclean_str = unclean_str.replace(ch, '')
    for ch in ['/', '\\', '|', ' ', ',', '.']: unclean_str = unclean_str.replace(ch, '_')
    for ch in ['=', '+', ':', ';']: unclean_str = unclean_str.replace(ch, '-')
    return unclean_str

#def get_kpts(atoms_obj):
#    structure = AseAtomAdaptor.get_structure(atoms_obj)
#    kpts_dat = get_kpoints(structure, minDistance = 30, include_gamma = False)
#    return kpts_dat['cords']


def secant_method(func: Callable[[float | int], float | int], guess_minus: float | int, guess_current: float | int, maxs_iter: int = 300, con_cri: float | int = 10**(-10)):
    func_eval_minus = func(guess_minus)
    nr_iter = 1

#    BASING THE SECANT METHOD ON THE FORWARD EULER OF THE FUNCTION RESULTS IN TO MUCH FLUCTUATION FOR CONVERGENCE SINCE THE STEPS OF THE SECANT METHOD ARE TO BROAD FOR THE FORWARD EULER TO A GOOD APPROXIMATION
#    forward_euler = lambda x_cur,x_minus,eval_cur,eval_minus: (eval_cur - eval_minus)/(x_cur-x_minus)
#    func_div_eval_minus = forward_euler(guess_minus, guess_minus * 0.95, func_eval_minus, func(guess_minus * 0.95))
#    while abs(func_div_eval_current := forward_euler(guess_current,guess_minus, func_eval_cur := func(guess_current),func_eval_minus)) >= con_cri and nr_iter < maxs_iter:
#        guess_current_temp = guess_current
#        guess_current -= (func_div_eval_current * (guess_current - guess_minus)) / (func_div_eval_current - func_div_eval_minus)
#        guess_minus = guess_current_temp
#        func_eval_minus = func_eval_cur
#        func_div_eval_minus = func_div_eval_current
#        nr_iter += 1
#    return guess_current, nr_iter

    while abs(func_eval_current := func(guess_current)) >= con_cri and nr_iter < maxs_iter:
        guess_current_temp = guess_current
        guess_current -= (func_eval_current * (guess_current - guess_minus)) / (func_eval_current - func_eval_minus)
        guess_minus = guess_current_temp
        func_eval_minus = func_eval_current
        nr_iter += 1
    return guess_current, nr_iter


def calculate_pE_of_latt_vasp(lattice: float, metal: str, slab_type: str, functional: str, functional_folder: str, grid_spacing: float, correction: Optional[str] = None, step_obj: Optional['steps'] = None) -> float:
    if (isinstance(lattice, list) or isinstance(lattice, tuple)) and len(lattice) == 1: lattice = lattice[0]
    lattice = float(lattice)


    match slab_type.lower():
        case 'fcc' | 'hcp' | 'bcc': atoms = bulk(name=metal, crystalstructure=slab_type, a=lattice)
        case 'graphene': atoms = graphene(metal, a=lattice, vacuum=10)
        case _: raise ValueError(f'Slab type was not recognised: {slab_type}')


#    match correction:
#        case 'DFTD4' | 'D4':
#            correction_func = DFTD4
#            correction_arg = dict(
#                method=functional
#            )
#        case _:
#            raise ValueError('Corrections were not recognised and have likely not been implement.')

    atoms.pbc = True

    calc = Vasp(atoms=atoms,
                txt=f'{functional_folder}/{metal}_latt_fit/lat-opt_{metal}_{slab_type}_a-{lattice}.txt',
                xc=functional,
                kpts=[11, 11, 11],
                system=f'{functional_folder}/{metal}_latt_fit/lat-opt_{metal}_{slab_type}_a-{lattice}.txt',  # Insert the name of the system
                istart=0,  # Wavefunction
                icharg=2,  # Charge: 1-file 2-atom 10-const
                # PREC = Normal           # Specifies the "precision"-mode
                #ispin=spin,  # Specifies spin polarization (1 no, 2 yes)
                encut=450,  # Energy limit for plane waves in eV (never <300)
                nelmin=6,  # Minimum number of electronic steps
                nelm=500,
                ismear=0,  # Smearing
                sigma=0.2,  # Smearing's Width (10^-3 for molecule)
                ibrion=2,  # Atomic optimization algorithm
                nsw=500,  # Maximum number of ionic steps
                potim=0.25,  # IBRION's scaling factor
                ediff=0.00001,  # Electronic optimization criteria (energy, eV)
                ediffg=-0.01,  # Atomic optimization criteria (negative: forces eV/A, positive: energy eV)
                # npar=4,  # number of bands that are treated in parallel, sqrt of n of cores
                lcharg=False,
                lwave=False,
                algo='A'
    )

    atoms.calc = calc

    potential_energy = atoms.get_potential_energy()
    ### DELETE BULK AND CALC ###
    del calc
    del atoms

    if world.rank == 0 and step_obj is not None:
        step_obj.lattice.append(lattice)
        step_obj.energy.append(potential_energy)

    return potential_energy



def report(res: OptimizeResult) -> None:
    parprint(f'optimisation: {"succesfull" if res.success else "unsuccesfull"}')
    parprint(f'finale lattice: {res.x}')


#def plot_steps(steps: 'steps', save_name: Optional[str]) -> None:
#    fig = go.Figure()
#    fig.add_trace(go.Scatter(
#                     x=steps.lattice, y=steps.energy,
#                     mode='markers+lines',
#                     line=dict(color='grey')
#    ))

#    fig.update_layout(
#        xaxis_title=f'lattice constant',
#        yaxis_title='potential energy eV'
#    )

#    if save_name: fig.write_html(save_name)#, include_mathjax='cdn')
#    else: fig.show()


def main(metal: str, functional: str, slab_type: str, guess_lattice: Optional[float] = None, grid_spacing: float = 0.16, correction: Optional[str] = None):

    functional_folder = sanitize('_'.join([functional, correction]) if correction is not None else functional)

    script_overlab_protection_time = randint(0, 60)
    if world.rank == 0:
        sleep(script_overlab_protection_time)
        folder_exist(functional_folder)
        folder_exist(f'{functional_folder}/{metal}_latt_fit')
    else:
        sleep(script_overlab_protection_time)

    if guess_lattice is None:
        if slab_type == 'graphene': guess_lattice = 2.46
        else:
            at_number = chemical_symbols.index(metal)
            if slab_type != reference_states[at_number].get('symmetry'): raise ValueError('the given slab type does not match the saved type for ase guess lattice')
            guess_lattice = reference_states[at_number].get('a')

    parprint(f'lattice optimisation for {metal} with {functional}, guess latice is at {guess_lattice}')


    opts_steps = make_dataclass('steps',
                           [('lattice', list, field(default_factory=list)),
                            ('energy', list, field(default_factory=list))]
    )()  # this obj is made to save the optimisation steps without having to return the values.

    opt_step_func = lambda lat: calculate_pE_of_latt_vasp(lat, metal, slab_type, functional, functional_folder, grid_spacing, correction=correction, step_obj=opts_steps) # this is to make a function which is only dependent on a single variable lat

#    optimised_lat,final_itr = secant_method(opt_step_func,guess_minus= guess_lattice*0.9, guess_current=guess_lattice,maxs_iter=30)
    opt_res = minimize(opt_step_func, x0=guess_lattice, method='Powell', tol=0.01, options=dict(disp=True), bounds=((1.5, 7),))

    report(opt_res)

    if world.rank == 0 and opt_res.success:
        if 'lattice_calc.csv' not in os.listdir(): pathlib.Path('lattice_calc.csv').touch()
        with open('lattice_calc.csv', 'a') as csv_file:
            fields = ['metal', 'type', 'functional', 'lattice']
            writer_obj = csv.DictWriter(csv_file, fieldnames=fields)
            writer_obj.writerow(
                dict(
                    metal=metal,
                    type=slab_type,
                    functional='_'.join([functional, correction]) if correction is not None else functional,
                    lattice=opt_res.x if not isinstance(opt_res.x, list) else opt_res.x[0]
                )
            )

    parprint(opts_steps)
#    if world.rank == 0: plot_steps(opts_steps, f'{functional_folder}/opt_steps_{metal}_{functional}.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('metal', type=str)
    parser.add_argument('surface_type', type=str, choices=('fcc', 'bcc', 'hcp', 'graphene'))
    parser.add_argument('func', type=str)
    parser.add_argument('--lattice', '-a', type=float)
    parser.add_argument('--correction', '-cor', type=str, choices=('DFTD4', 'D4'))
    args = parser.parse_args()

    main(metal=args.metal, functional=args.func, slab_type=args.surface_type, guess_lattice=args.lattice, correction=args.correction)

