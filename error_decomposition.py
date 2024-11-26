import sys
import pathlib
#from copy import copy
from itertools import chain

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
#from error_project_san_sebastion import build_pd
#from error_project_san_sebastion.reaction_functions import Functional, reaction, get_needed_structures
from reaction_functions import Functional, reaction, get_needed_structures
#from error_project_san_sebastion.reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named
from reactions import all_gaseous_reactions, all_formation_reactions, all_gaseous_reactions_named, all_formation_reactions_named

import numpy as np
import pandas as pd
from scipy.linalg import lstsq


def build_error_matrix(reactions: list[reaction], functional_groups: dict[str, dict[str, int]]) -> pd.DataFrame:
    all_components = {comp.name for comp in chain.from_iterable(comp_list for comp_list in chain.from_iterable((react.reactants, react.products) for react in reactions))}
    all_functional_groups = tuple({str(group) for group in chain.from_iterable(functional_groups[comp].keys() for comp in all_components if comp in functional_groups.keys())})

    reactions_group_matrix = pd.DataFrame(dict(
        reaction=reactions,
        **{group: [0]*len(reactions) for group in all_functional_groups}
        ))

    def set_group_values(seq: pd.Series):
        seq_copy = seq.copy()
        for reactant in seq['reaction'].reactants:
            if reactant.name in functional_groups.keys():
                for group, number in functional_groups[reactant.name].items():
                    seq_copy[group] -= number * reactant.amount
        for product in seq['reaction'].products:
            if product.name in functional_groups.keys():
                for group, number in functional_groups[product.name].items():
                    seq_copy[group] += number * product.amount
        return seq_copy

    reactions_group_matrix = reactions_group_matrix.apply(set_group_values, axis=1)
    return reactions_group_matrix


def lstsq_decomposition(functional_obj: Functional, reactions: list[reaction], functional_groups: dict[str, dict[str, int]]) -> tuple[dict[str, float], float]:
    deviations = np.array([None] * len(reactions))
    for i, reac in enumerate(reactions):
        try: deviations[i] = functional_obj.calculate_reaction_enthalpy(reac) - reac.experimental_ref
        except: pass

    decomposition_matrix = build_error_matrix(reactions, functional_groups)
    decomposition_matrix, deviations = remove_nan_error_matrix(decomposition_matrix, deviations)

    group_errors, residual, _, _ = lstsq(decomposition_matrix.drop('reaction', axis=1), deviations)
    errors = {col: err for col, err in zip(decomposition_matrix.drop('reaction', axis=1).columns, group_errors)}

    return errors, residual


def remove_nan_error_matrix(decomposition_matrix: pd.DataFrame, deviation_vector: np.typing.ArrayLike) -> tuple[pd.DataFrame, np.typing.ArrayLike]:
    decomposition_matrix_filtered = decomposition_matrix.loc[~np.isnan(deviation_vector)] #lambda x: deviation_vector[x['id']] is not None)
    deviation_vector_filtered = deviation_vector[~np.isnan(deviation_vector)]
    return decomposition_matrix_filtered, deviation_vector_filtered


def simple_decomposition(functional_obj: Functional, functional_groups: dict[str, dict[str, int]]) -> dict[str, float]:
    reaction_groups = (
        'Simple alkanes',
        'Branched alkanes (iso)',
        'Branched alkanes (neo)',
        'Amines',
        'Nitros',
        'Nitrates',
        'Nitrites',
        'Hydroxylamines',
        'Aromatics',
        'Anilines',
        'Hydrazines',
        'Amides',
        'Nitriles',
    )

    error_foci = (
        'methyl_carbons',
        'iso_carbons',
        'neo_carbons',
        'amines',
        'nitro',
        'nitrate',
        'nitrite',
        'hydroxylamine',
        'phenyl',
        'aniline',
        'hydrazine',
        'amide',
        'nitrile',
    )

    errors = {}
    for err_foc, chem_name in zip(error_foci, reaction_groups):
        errors.update({
            err_foc: np.mean([(functional_obj.calculate_reaction_enthalpy(reac) - reac.experimental_ref - sum(errors[err] * functional_groups[reac.product[0]][err] for err in errors.keys() if err in functional_groups.keys())) / functional_groups[reac.product[0]][err_foc] for reac in all_gaseous_reactions_named[chem_name]])
        })
    return errors
