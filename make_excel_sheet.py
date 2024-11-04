import argparse
import sys
import pathlib
from copy import copy

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import build_pd
from error_project_san_sebastion.reaction_functions import Functional, get_needed_structures
from error_project_san_sebastion.reactions import all_gaseous_reactions, all_formation_reactions

import pandas as pd
import openpyxl as xl
from openpyxl.formatting.rule import ColorScaleRule


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



def main(molecule_database_dir: str, solid_database_dir: str, verbose: bool = False):
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

    excel_file = xl.Workbook()
    work_sheet = excel_file.active
    work_sheet.title = 'energies'
    deviation_sheet = excel_file.create_sheet('deviations')
    formation_sheet = excel_file.create_sheet('formation')
    formation_sheet_deviation = excel_file.create_sheet('formation_deviation')
    correction_sheet = excel_file.create_sheet('corrections')

    for i, func in enumerate(functional_list):
        for sheet in (work_sheet, deviation_sheet, formation_sheet, formation_sheet_deviation, correction_sheet):
            sheet.cell(1, i+2, func.name)
            correction_sheet.cell(2, i+2, O2_er[func.name])
            correction_sheet.cell(3, i+2, N2_er[func.name])
    correction_sheet.cell(2, 1, 'O2 error')
    correction_sheet.cell(3, 1, 'N2 error')

    work_sheet.cell(1, len(functional_list) + 2, 'exp ref')
    formation_sheet.cell(1, len(functional_list) + 2, 'exp ref')

    for i, reac in enumerate(all_gaseous_reactions):
        work_sheet.cell(i+2, 1, reac.products[0].name)
        deviation_sheet.cell(i+2, 1, reac.products[0].name)
        for j, func in enumerate(functional_list):
            try:
                work_sheet.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac))
                deviation_sheet.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac) - reac.experimental_ref)
            except: pass
        work_sheet.cell(i+2, len(functional_list) + 2, reac.experimental_ref)
    deviation_sheet.conditional_formatting.add(f'B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}',
                                               ColorScaleRule(start_type='formula',
                                                              start_value=f'=-MAX(-MIN(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2});MAX(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}))',
                                                              start_color='AA0000',

                                                              mid_type='num',
                                                              mid_value=0,
                                                              mid_color='FFFFFF',

                                                              end_type='formula',
                                                              end_value=f'=-MAX(-MIN(B2:{xl.utils.cell.get_column_letter(j)}{i});MAX(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}))',
                                                              end_color='AA0000'
                                                              ))

    for i, reac in enumerate(all_formation_reactions):
        formation_sheet.cell(i+2, 1, reac.products[0].name)
        formation_sheet_deviation.cell(i+2, 1, reac.products[0].name)
        for j, func in enumerate(functional_list):
            try:
                formation_sheet.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac))
                formation_sheet_deviation.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac) - reac.experimental_ref)
            except: pass
        formation_sheet.cell(i+2, len(functional_list) + 2, reac.experimental_ref)
    formation_sheet_deviation.conditional_formatting.add(f'B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}',
                                               ColorScaleRule(start_type='formula',
                                                              start_value=f'=-MAX(-MIN(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2});MAX(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}))',
                                                              start_color='AA0000',

                                                              mid_type='num',
                                                              mid_value=0,
                                                              mid_color='FFFFFF',

                                                              end_type='formula',
                                                              end_value=f'=-MAX(-MIN(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2});MAX(B2:{xl.utils.cell.get_column_letter(j+2)}{i+2}))',
                                                              end_color='AA0000'
                                                              ))

    excel_file.save('reaction_results.xlsx')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('molecule_db_directory', help='Path to the database')
    parser.add_argument('solid_db_directory', help='Path to the database')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    main(args.molecule_db_directory, args.solid_db_directory, args.verbose)

