import argparse
import sys
import pathlib
from copy import copy

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
#from error_project_san_sebastion import build_pd
from . import build_pd
#from error_project_san_sebastion.reaction_functions import Functional, get_needed_structures
from reaction_functions import Functional, get_needed_structures
#from error_project_san_sebastion.reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named
from reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named
#from error_project_san_sebastion.error_decomposition import simple_decomposition
from error_decomposition import simple_decomposition, lstsq_decomposition
#from error_project_san_sebastion.manual_functional_groups import molecule_functional_dict
from manual_functional_groups_revised import molecule_functional_dict

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

    upper_border = Border(
        top=Side(border_style=BORDER_THIN, color='00000000'),
        )

    start_of_data = 3

    for i, func in enumerate(functional_list):
        for sheet in (work_sheet, deviation_sheet, formation_sheet, formation_sheet_deviation): sheet.cell(1, i + start_of_data, func.name)
        correction_sheet.cell(1, i + start_of_data - 1, func.name)
        correction_sheet.cell(2, i+2, O2_er[func.name])
        correction_sheet.cell(3, i+2, N2_er[func.name])
    correction_sheet.cell(2, 1, 'O2 error')
    correction_sheet.cell(3, 1, 'N2 error')

    work_sheet.cell(1, len(functional_list) + start_of_data, 'exp ref')
    formation_sheet.cell(1, len(functional_list) + start_of_data, 'exp ref')

    next_row = 2

    for name, reactions in all_gaseous_reactions_named.items():
        for sheet in [work_sheet, deviation_sheet]:
            sheet.cell(next_row, start_of_data - 2, name)
            for j in range(len(functional_list)+start_of_data+1): sheet.cell(next_row, j + 1).border = upper_border
        for i, reac in enumerate(reactions):
            for sheet in [work_sheet, deviation_sheet]: sheet.cell(next_row, start_of_data - 1, reac.products[0].name)
            for j, func in enumerate(functional_list):
                try:
                    work_sheet.cell(next_row, j + start_of_data, func.calculate_reaction_enthalpy(reac))
                    deviation_sheet.cell(next_row, j + start_of_data, func.calculate_reaction_enthalpy(reac) - reac.experimental_ref)
                except: pass
            work_sheet.cell(next_row, len(functional_list) + start_of_data, reac.experimental_ref)
            next_row += 1

    deviation_sheet.conditional_formatting.add((cond_zone := f'${xl.utils.cell.get_column_letter(start_of_data)}${start_of_data}:${xl.utils.cell.get_column_letter(j+start_of_data)}${next_row - 1}'),
                                               ColorScaleRule(start_type='formula',
                                                              start_value=f'-MAX(-MIN({cond_zone}),MAX({cond_zone}))',
                                                              start_color='AA0000',

                                                              mid_type='num',
                                                              mid_value=0,
                                                              mid_color='FFFFFF',

                                                              end_type='formula',
                                                              end_value=f'MAX(-MIN({cond_zone}),MAX({cond_zone}))',
                                                              end_color='AA0000'
                                                              ))

    next_row = 2

    for name, reactions in all_formation_reactions_named.items():
        for sheet in [formation_sheet, formation_sheet_deviation]:
            sheet.cell(next_row, start_of_data - 2, name)
            for j in range(len(functional_list)+start_of_data+1): sheet.cell(next_row, j + 1).border = upper_border
        for i, reac in enumerate(reactions):
            for sheet in [formation_sheet, formation_sheet_deviation]: sheet.cell(next_row, start_of_data - 1, reac.products[0].name)
            for j, func in enumerate(functional_list):
                try:
                    formation_sheet.cell(next_row, j + start_of_data, func.calculate_reaction_enthalpy(reac))
                    formation_sheet_deviation.cell(next_row, j + start_of_data, func.calculate_reaction_enthalpy(reac) - reac.experimental_ref)
                except: pass
            formation_sheet.cell(next_row, len(functional_list) + start_of_data, reac.experimental_ref)
            next_row += 1
    formation_sheet_deviation.conditional_formatting.add((cond_zone := f'${xl.utils.cell.get_column_letter(start_of_data)}${start_of_data}:${xl.utils.cell.get_column_letter(j+start_of_data)}${next_row - 1}'),
                                               ColorScaleRule(start_type='formula',
                                                              start_value=f'-MAX(-MIN({cond_zone}),MAX({cond_zone}))',
                                                              start_color='AA0000',

                                                              mid_type='num',
                                                              mid_value=0,
                                                              mid_color='FFFFFF',

                                                              end_type='formula',
                                                              end_value=f'MAX(-MIN({cond_zone}),MAX({cond_zone}))',
                                                              end_color='AA0000'
                                                              ))

    reac_group_names = {'Simple alkanes': 'methyl_carbons', 'Branched alkanes (iso)': 'iso_carbons', 'Branched alkanes (neo)': 'neo_carbons', 'Amines': 'amines', 'Nitros': 'nitro', 'Nitrates': 'nitrate', 'Nitrites': 'nitrite', 'Hydroxylamines': 'hydroxylamine', 'Aromatics': 'phenyl', 'Anilines': 'aniline', 'Hydrazines': 'hydrazine', 'Amides': 'amide', 'Nitriles': 'nitrile'}

    for i, func in enumerate(functional_list):
        try:
            corrections = simple_decomposition(func, molecule_functional_dict)
            for j, (group_name, val) in enumerate(corrections.items()):
                if correction_sheet.cell(1, 5 + j).value is None: correction_sheet.cell(1, 5 + j, group_name)
                correction_sheet.cell(2 + i, 5 + j, val)
        except: pass

    correction_sheet_linalg = excel_file.create_sheet('corrections linalg')
    correction_sheet_linalg.cell(1, 2, 'residual')
    start_of_linalg_data = 5

    for i, func in enumerate(functional_list):
        try:
            correction_sheet_linalg.cell(1 + i, 1, func.name)
            corrections, residual = lstsq_decomposition(func, all_gaseous_reactions, molecule_functional_dict)
            correction_sheet_linalg.cell(1 + i, 2, residual)
            for j, (group_name, corr_key) in enumerate(reac_group_names.items()):
                if correction_sheet_linalg.cell(1, start_of_linalg_data + j).value is None: correction_sheet_linalg.cell(1, start_of_linalg_data + j, group_name)
                correction_sheet_linalg.cell(2 + i, start_of_linalg_data + j, corrections[corr_key])
        except: pass


    excel_file.save('reaction_results.xlsx')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('molecule_db_directory', help='Path to the database')
    parser.add_argument('solid_db_directory', help='Path to the database')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    main(args.molecule_db_directory, args.solid_db_directory, args.verbose)
