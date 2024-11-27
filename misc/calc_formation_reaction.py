import argparse
import sys
import pathlib
from copy import copy
import traceback

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.parent))
from error_project_san_sebastion import build_pd
#from . import build_pd
from error_project_san_sebastion.reaction_functions import Functional, get_needed_structures
#from .reaction_functions import Functional, get_needed_structures
from error_project_san_sebastion.reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named
#from .reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named
from error_project_san_sebastion.error_decomposition import simple_decomposition, lstsq_decomposition
#from .error_decomposition import simple_decomposition, lstsq_decomposition
from error_project_san_sebastion.manual_functional_groups_revised import molecule_functional_dict
#from .manual_functional_groups_revised import molecule_functional_dict

import pandas as pd
import openpyxl as xl
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles.borders import Border, Side, BORDER_THIN


def missing_structures(functional: Functional, needed_strcutures: dict):
    missing_molecules = {key for key in needed_strcutures['molecule'] if key not in functional.molecule.keys()}
    missing_adsorbates = {key for key in needed_strcutures['adsorbate'] if key not in functional.adsorbate.keys()}
    missing_slab = {key for key in needed_strcutures['slab'] if key not in functional.slab.keys()}

    if any(len(mis)>0 for mis in (missing_adsorbates, missing_slab, missing_molecules)):
        print(f'{functional.name} is missing')
        if len(missing_molecules) > 0:
            print('    molecules:')
            print(missing_molecules)
        if len(missing_adsorbates) > 0:
            print('    adsorbates:')
            print(missing_adsorbates)
        if len(missing_slab) > 0:
            print('    slab:')
            print(missing_slab)
        print()


def main(molecule_database_dir: str, solid_database_dir: str, target_molecule: str, verbose: bool = False):
    pd_molecule_dat = build_pd(molecule_database_dir)
    pd_solid_dat = build_pd(solid_database_dir)

    pd_molecule_dat['enthalpy'] = pd_molecule_dat['energy'] + pd_molecule_dat['zpe']
    pd_solid_dat['enthalpy'] = pd_solid_dat['energy'] + pd_solid_dat['zpe']

    functional_set = {xc for _, row in pd_molecule_dat.iterrows() if not pd.isna((xc := row.get('xc')))}

    need_structures = get_needed_structures(all_gaseous_reactions + all_formation_reactions)

    functional_list = []
    for xc in functional_set:
        try: functional_list.append(Functional(functional_name=xc, slab_db=pd_solid_dat, adsorbate_db=None, mol_db=pd_molecule_dat, needed_struc_dict=need_structures, thermo_dynamic=True))
        except: pass

    if verbose:
        for i, func in enumerate(functional_list):
            missing_structures(func, need_structures)

    O2_er = dict()
    N2_er = dict()

    for func in copy(functional_list):
        try:
            O2_er.update({func.name: func.calc_O2_err()})
            N2_er.update({func.name: func.calc_N2_err()})
            func.correct_energies()
        except: functional_list.remove(func)

    react = list(filter(lambda r: r.products[0].name == target_molecule, all_formation_reactions))[0]
    result = {}
    reactants = {}
    products = {}

    for xc in functional_list:
        try:
            result.update({xc.name: xc.calculate_reaction_enthalpy(react)})
            reactants.update({rea.name: getattr(xc, rea.type)[rea.name] for rea in react.reactant})
            products.update({pro.name: getattr(xc, pro.type)[pro.name] for pro in react.reactant})
        except: pass

    print(f'formation energy is: {result}')
    print(f'reactant energies are {reactants}')
    print(f'product energies are {products}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('molecule_db_directory', help='Path to the database')
    parser.add_argument('solid_db_directory', help='Path to the database')
    parser.add_argument('molecule', help='Path to the database')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    main(args.molecule_db_directory, args.solid_db_directory, args.molecule, args.verbose)
