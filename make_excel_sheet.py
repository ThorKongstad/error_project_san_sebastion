import argparse
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion import build_pd
from error_project_san_sebastion.reaction_functions import Functional, get_needed_structures
from error_project_san_sebastion.reactions import all_reactions

import pandas as pd
import openpyxl as xl


def main(database_dir: str):
    pd_dat = build_pd(database_dir)

    functional_set = {xc for _, row in pd_dat.iterrows() if not pd.isna((xc := row.get('xc')))}

    need_structures = get_needed_structures(all_reactions)

    functional_list = []
    for xc in functional_set:
        try: functional_list.append(Functional(functional_name=xc, slab_db=None, adsorbate_db=None, mol_db=pd_dat, needed_struc_dict=need_structures, thermo_dynamic=False))
        except: pass

    excel_file = xl.Workbook()
    work_sheet = excel_file.active
    work_sheet.title = 'energies'
    deviation_sheet = excel_file.create_sheet('deviations')

    for i, func in enumerate(functional_list):
        work_sheet.cell(1, i+2, func.name)
        deviation_sheet.cell(1, i+2, func.name)
    work_sheet.cell(1, len(functional_list) + 2, 'exp ref')

    for i, reac in enumerate(all_reactions):
        work_sheet.cell(i+2, 1, reac.products[0].name)
        deviation_sheet.cell(i+2, 1, reac.products[0].name)
        for j, func in enumerate(functional_list):
            work_sheet.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac))
            deviation_sheet.cell(i+2, j+2, func.calculate_reaction_enthalpy(reac) - reac.experimental_ref)
        work_sheet.cell(i+2, len(functional_list) + 2, reac.experimental_ref)

    excel_file.save('reaction_results.xlsl')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('db_directory', help='Path to the database')
    args = parser.parse_args()

    main(args.db_directory)

