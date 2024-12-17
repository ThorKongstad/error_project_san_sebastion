from typing import NoReturn, Sequence, Tuple, Never, Optional, NamedTuple, Iterable
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from ase.db.core import bytes_to_object



class component(NamedTuple):
    type: str
    name: str
    amount: float


@dataclass
class reaction:
    reactants: Tuple[Tuple[str, str, float] | component, ...]
    products: Tuple[Tuple[str, str, float] | component, ...]
    experimental_ref: float | None = None

    def __post_init__(self):
        for n, reacs_or_prods in enumerate([self.reactants, self.products]):
            new_component_seq = []
            for i, reac_or_prod in enumerate(reacs_or_prods):
                if len(reac_or_prod) != 3: raise ValueError('a component of a reaction does not have the correct size')
                if not reac_or_prod[0] in ('molecule', 'slab', 'adsorbate'): raise ValueError('The reactant or product type string appear to be wrong')
                new_component_seq.append(component(*reac_or_prod) if not isinstance(reac_or_prod, component) else reac_or_prod)
            setattr(self, 'reactants' if n == 0 else 'products', tuple(new_component_seq))

    def __str__(self):
        return ' ---> '.join([' + '.join([f'{reac.amount:.2g}{reac.name if reac.name != "cid281" else "C|||O"}({reac.type})' for reac in comp]) for comp in (self.reactants, self.products)])


class Functional:
    def __init__(self, functional_name: str, slab_db: pd.DataFrame, adsorbate_db: pd.DataFrame, mol_db: pd.DataFrame, needed_struc_dict: Optional[dict[str, Iterable[str]]] = None, thermo_dynamic: bool = True):
        energy_type = 'enthalpy' if thermo_dynamic else 'energy'
        self.name = functional_name

        self.molecule = {}
        if mol_db is not None:
            for smile in needed_struc_dict['molecule']:
                try: self.molecule.update({smile: mol_db.query(f'name == "{smile}" and xc == "{functional_name}" and {energy_type}.notna()').get(energy_type).iloc[0]})
                except: pass

        self.slab = {}
        if slab_db is not None:
            for structure_str in needed_struc_dict['slab']:
                try: self.slab.update({structure_str: slab_db.query(f'name == "{structure_str}" and xc == "{functional_name}" and {energy_type}.notna()').get(energy_type).iloc[0]})
                except: pass

        self.adsorbate = {}
        if adsorbate_db is not None:
            for structure_str in needed_struc_dict['adsorbate']:
                try: self.adsorbate.update({structure_str: adsorbate_db.query(f'structure_str == "{structure_str}" and xc == "{functional_name}" and {energy_type}.notna()').get(energy_type).iloc[0]})
                except: pass

        self.has_BEE = functional_name == 'BEEF-vdW'
        if self.has_BEE:

            self.molecule_energy = {}
            self.molecule_bee = {}
            if mol_db is not None:
                for smile in needed_struc_dict['molecule']:
                    try: self.molecule_energy.update({smile: mol_db.query(f'name == "{smile}" and xc == "{functional_name}" and energy.notna()').get('energy').iloc[0]})
                    except: pass
                    try: self.molecule_bee.update({smile: np.array(bytes_to_object(mol_db.query(f'name == "{smile}" and xc == "{functional_name}" and _data.notna()').get('_data').iloc[0]).get('ensemble_en'))[:]})
                    except: pass


            self.slab_energy = self.slab
            self.slab_bee = {}
            for structure_str in needed_struc_dict['slab']:
                try: self.slab_bee.update({structure_str: np.array(bytes_to_object(slab_db.query(f'name == "{structure_str}" and xc == "{functional_name}" and _data.notna()').get('_data').iloc[0]).get('ensemble_en'))[:]})
                except: pass


            self.adsorbate_energy = {}
            self.adsorbate_bee = {}
            if adsorbate_db is not None:
                for structure_str in needed_struc_dict['adsorbate']:
                    try: self.adsorbate_energy.update({structure_str: adsorbate_db.query(f'structure_str == "{structure_str}" and xc == "{functional_name}" and energy.notna()').get('energy').iloc[0]})
                    except: pass
                    try: self.adsorbate_bee.update({structure_str: np.array(bytes_to_object(adsorbate_db.query(f'structure_str == "{structure_str}" and xc == "{functional_name}" and _data.notna()').get('_data').iloc[0]).get('ensemble_en'))[:]})
                    except: pass

        else:
            self.molecule_energy = {}
            self.adsorbate_energy ={}
            self.slab_energy = {}
            self.molecule_bee = {}
            self.slab_bee = {}
            self.adsorbate_bee = {}

    def calculate_reaction_enthalpy(self, reaction: reaction) -> float:
        reactant_enthalpy, product_enthalpy = tuple(sum(getattr(self, typ)[name] * amount for typ, name, amount in getattr(reaction, reac_part)) for reac_part in ('reactants', 'products'))
        return product_enthalpy - reactant_enthalpy

    def calculate_reaction_energy(self, reaction: reaction) -> float:
        if not self.has_BEE: raise ValueError('calculate_reaction_energy only if the functional has BEE')
        reactant_enthalpy, product_enthalpy = tuple(sum(getattr(self, typ + '_energy')[name] * amount for typ, name, amount in getattr(reaction, reac_part)) for reac_part in ('reactants', 'products'))
        return product_enthalpy - reactant_enthalpy

    def calculate_BEE_reaction_enthalpy(self, reaction: reaction) -> np.ndarray[float]:
        if not self.has_BEE: raise ValueError('calculate_reaction_energy only if the functional has BEE')
        correction = self.calculate_reaction_enthalpy(reaction) - self.calculate_reaction_energy(reaction)
        reactant_BEE_enthalpy, product_BEE_enthalpy = tuple(sum(getattr(self, typ + '_bee')[name] * amount for typ, name, amount in getattr(reaction, reac_part)) for reac_part in ('reactants', 'products'))
        return product_BEE_enthalpy - reactant_BEE_enthalpy + correction

    def calc_N2_err(self) -> float:
        amonium_reaction = reaction((('molecule', 'nitrogen', 0.5), ('molecule', 'hydrogen', 1.5),), (('molecule', 'ammonia', 1),), -0.48)

        return -2 * (self.calculate_reaction_enthalpy(amonium_reaction) - amonium_reaction.experimental_ref)

    def calc_O2_err(self) -> float:
#        peroxide_reaction = reaction((('molecule', 'oxygen', 1), ('molecule', 'hydrogen', 1),), (('molecule', 'peroxide', 1),), -1.09)
        water_reaction = reaction((('molecule', 'oxygen', 0.5), ('molecule', 'hydrogen', 1),), (('molecule', 'water', 1),), -2.51)

        return -2 * (self.calculate_reaction_enthalpy(water_reaction) - water_reaction.experimental_ref)

    def calc_CO_err(self) -> float:
    #        peroxide_reaction = reaction((('molecule', 'oxygen', 1), ('molecule', 'hydrogen', 1),), (('molecule', 'peroxide', 1),), -1.09)
        co_reaction = reaction((('molecule', 'oxygen', 0.5), ('slab', 'graphene', 0.5),), (('molecule', 'carbon-monoxide', 1),), -1.15)

        return (self.calculate_reaction_enthalpy(co_reaction) - co_reaction.experimental_ref)

    def calc_NO_err(self) -> float:
    #        peroxide_reaction = reaction((('molecule', 'oxygen', 1), ('molecule', 'hydrogen', 1),), (('molecule', 'peroxide', 1),), -1.09)
        co_reaction = reaction((('molecule', 'oxygen', 0.5), ('molecule', 'nitrogen', 0.5),), (('molecule', 'nitric-oxide', 1),), 0.95)

        return (self.calculate_reaction_enthalpy(co_reaction) - co_reaction.experimental_ref)


    def correct_energies(self):
        self.molecule['oxygen'] -= self.calc_O2_err()
        self.molecule['nitrogen'] -= self.calc_N2_err()
        self.molecule['nitric-oxide'] -= self.calc_NO_err()
        self.molecule['carbon-monoxide'] -= self.calc_CO_err()


def get_needed_structures(reactions_seq: Sequence[reaction]):
    dictionary_of_needed_strucs = {'molecule': [], 'slab': [], 'adsorbate': []}
    for reac in reactions_seq:
        for compo in reac.reactants + reac.products:
            if compo.name not in dictionary_of_needed_strucs[compo.type]: dictionary_of_needed_strucs[compo.type].append(compo.name)
    return dictionary_of_needed_strucs