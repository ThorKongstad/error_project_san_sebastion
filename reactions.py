import sys
import pathlib
from typing import Tuple

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion.reaction_functions import reaction

#reaction((('',),), (('',),), experimental value)
# https://doi.org/10.1002/cctc.202100125


simple_alkane_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 2),), (('molecule', 'methane', 1),), -0.77),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 3),), (('molecule', 'ethane', 1),), -0.87),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 4),), (('molecule', 'propane', 1),), -1.08),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5),), (('molecule', 'butane', 1),), -1.30),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 6),), (('molecule', 'pentane', 1),), -1.52),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7),), (('molecule', 'hexane', 1),), -1.73),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 8),), (('molecule', 'heptane', 1),), -1.94),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9),), (('molecule', 'octane', 1),), -2.16),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 11),), (('molecule', 'decane', 1),), -2.59),
)

iso_alkane_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5),), (('molecule', 'isobutane', 1),), -1.39),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 6),), (('molecule', 'isopentane', 1),), -1.59),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7),), (('molecule', 'isohexane', 1),), -1.81),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7),), (('molecule', '3-methylpentane', 1),), -1.78),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 8),), (('molecule', 'isoheptane', 1),), -2.02),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 8),), (('molecule', '3-methylhexane', 1),), -1.98),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9),), (('molecule', '2-methylheptane', 1),), -2.23),
)

neo_alakane_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 6),), (('molecule', 'neopentane', 1),), -1.74),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7),), (('molecule', 'neohexane', 1),), -1.93),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 8),), (('molecule', 'neoheptane', 1),), -2.13),
)

amine_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitrogen', 0.5), ('molecule', 'hydrogen', 1.5),), (('molecule', 'ammonia', 1),), -0.48),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'methylamine', 1),), -0.23),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'ethylamine', 1),), -0.49),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'propylamine', 1),), -0.73),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'isopropylamine', 1),), -0.87),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', '1-butylamine', 1),), -0.95),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', '2-butylamine', 1),), -1.08),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'isobutylamine', 1),), -1.02),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'tert-butylamine', 1),), -1.25),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'dimethylamine', 1),), -0.19),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'diethylamine', 1),), -0.75),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 6.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'methyltertbutylamine', 1),), -1.12),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'dipropylamine', 1),), -1.2),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'diisopropylamine', 1),), -1.49),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 8.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'isopropyl-tertbutylamine', 1),), -1.71),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'dibutylamine', 1),), -1.62),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'butylisobutylamine', 1),), -1.77),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'diisobutylamine', 1),), -1.86),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 9.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'ditertbutylamine', 1),), -1.78),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 7.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'triethylamine', 1),), -0.96),
    reaction((('slab', 'graphene', 4.5), ('molecule', 'hydrogen', 10.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'tripropylamine', 1),), -1.67),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'trimethylamine', 1),), -0.24),
)

nitro_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'nitromethane', 1),), -0.84),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'nitroethane', 1),), -1.08),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '1-nitropropane', 1),), -1.29),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '2-nitropropane', 1),), -1.44),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '1-nitrobutane', 1),), -1.49),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '2-nitrobutane', 1),), -1.7),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '1-nitropentane', 1),), -1.7),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', '2-methyl-2-nitropropane', 1),), -1.84),
)

nitrate_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'hydrogen', 0.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1.5)), (('molecule', 'nitric-acid', 1),), -1.39),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1.5)), (('molecule', 'methylnitrate', 1),), -1.26),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1.5)), (('molecule', 'ethylnitrate', 1),), -1.60),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1.5)), (('molecule', 'propylnitrate', 1),), -1.80),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1.5)), (('molecule', 'isopropylnitrate', 1),), -1.98),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 1.5), ('molecule', 'oxygen', 4.5)), (('molecule', '1.2.3-propanetriol-trinitrate', 1),), -2.89),
)


nitrite_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'hydrogen', 0.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'nitrous-acid-trans', 1),), -0.82),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'methyl-nitrite', 1),), -0.69),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'ethyl-nitrite', 1),), -1.12),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'propyl-nitrite', 1),), -1.23),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'isopropyl-nitrite', 1),), -1.38),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'butyl-nitrite', 1),), -1.51),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'sec-butyl-nitrite', 1),), -1.59),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 1)), (('molecule', 'tert-butylnitrite', 1),), -1.78),
)

hydroxylamine_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'hydroxylamine', 1),), -0.51),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'n.n-diethylhydroxylamine', 1),), -1.26),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'o-methylhydroxylamine', 1),), -0.26),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'n-methylhydroxylamine', 1),), -0.52),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'n.o-dimethylhydroxylamine', 1),), -0.39),
)

aromatic_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 3),), (('molecule', 'benzene', 1),), 0.86),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 4),), (('molecule', 'toluene', 1),), 0.52),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 5),), (('molecule', 'ethylbenzene', 1),), 0.31),
    reaction((('slab', 'graphene', 4.5), ('molecule', 'hydrogen', 6),), (('molecule', 'propylbenzene', 1),), 0.08),
    reaction((('slab', 'graphene', 4.5), ('molecule', 'hydrogen', 6),), (('molecule', 'isopropylbenzene', 1),), 0.04),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 7),), (('molecule', 'butylbenzene', 1),), -0.12),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 7),), (('molecule', 'sec-butylbenzene', 1),), -0.19),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 7),), (('molecule', 'tert-butylbenzene', 1),), -0.24),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 7),), (('molecule', 'isobutylbenzene', 1),), -0.23),
)

aniline_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'aniline', 1),), 0.91),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', '2-methylaniline', 1),), 0.58),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'n-methylaniline', 1),), 0.87),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'n-ethylaniline', 1),), 0.58),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'n.n-dimethylaniline', 1),), 1.04),
)

hydrazine_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'hydrogen', 2), ('molecule', 'nitrogen', 1)), (('molecule', 'hydrazine', 1),), 2.10),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 4), ('molecule', 'nitrogen', 1)), (('molecule', 'phenylhydrazine', 1),), 0.96),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 4), ('molecule', 'nitrogen', 1)), (('molecule', '1.2-dimethylhydrazine', 1),), 0.87),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 4), ('molecule', 'nitrogen', 1)), (('molecule', '1.1-dimethylhydrazine', 1),), 0.99),
)

amide_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'formamide', 1),), -2.01),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'acetamide', 1),), -2.47),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'propanamide', 1),), -2.68),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'butanamide', 1),), -2.92),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'pentanamide', 1),), -3.01),
    reaction((('slab', 'graphene', 3), ('molecule', 'hydrogen', 6.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'hexanamide', 1),), -3.36),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 8.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'octanamide', 1),), -3.76),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', 'isobutiramide', 1),), -2.93),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 5.5), ('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5)), (('molecule', '2.2-dimethylpropanamide', 1),), -3.25),
)

nitrile_formation: Tuple[reaction, ...] = (
    reaction((('slab', 'graphene', 0.5), ('molecule', 'hydrogen', 0.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'hydrogen-cyanide', 1),), 1.40),
    reaction((('slab', 'graphene', 1), ('molecule', 'hydrogen', 1.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'acetonitrile', 1),), 0.77),
    reaction((('slab', 'graphene', 1.5), ('molecule', 'hydrogen', 2.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'propanenitrile', 1),), 0.54),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'butanenitrile', 1),), 0.35),
    reaction((('slab', 'graphene', 2.5), ('molecule', 'hydrogen', 4.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'pentanenitrile', 1),), 0.11),
    reaction((('slab', 'graphene', 3.5), ('molecule', 'hydrogen', 6.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'heptanenitrile', 1),), -0.32),
    reaction((('slab', 'graphene', 4), ('molecule', 'hydrogen', 7.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'octanenitrile', 1),), -0.52),
    reaction((('slab', 'graphene', 5), ('molecule', 'hydrogen', 9.5), ('molecule', 'nitrogen', 0.5)), (('molecule', 'decanenitrile', 1),), -0.95),
    reaction((('slab', 'graphene', 2), ('molecule', 'hydrogen', 3.5), ('molecule', 'nitrogen', 0.5)), (('molecule', '2-methylpropanenitrile', 1),), 0.24),
)

small_molecule_formation: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitrogen', 0.5), ('molecule', 'oxygen', 0.5),), (('molecule', 'nitric-oxide', 1),), 0.95),
    reaction((('slab', 'graphene', 0.5), ('molecule', 'oxygen', 0.5),), (('molecule', 'carbon-monoxide', 1),), -1.15),
)

all_formation_reactions = simple_alkane_formation + iso_alkane_formation + neo_alakane_formation + amine_formation + nitro_formation + nitrate_formation + nitrite_formation + hydroxylamine_formation + aromatic_formation + aniline_formation + hydrazine_formation + amide_formation + nitrile_formation + small_molecule_formation
all_formation_reactions_named = {
    'Simple alkanes': simple_alkane_formation,
    'Branched alkanes (iso)': iso_alkane_formation,
    'Branched alkanes (neo)': neo_alakane_formation,
    'Amines': amine_formation,
    'Nitros': nitro_formation,
    'Nitrates': nitrate_formation,
    'Nitrites': nitrite_formation,
    'Hydroxylamines': hydroxylamine_formation,
    'Aromatics': aromatic_formation,
    'Anilines': aniline_formation,
    'Hydrazines': hydrazine_formation,
    'Amides': amide_formation,
    'Nitriles': nitrile_formation,
    'Small molecules': small_molecule_formation
}

####
simple_alkane_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'water', 2)), (('molecule', 'methane', 1), ('molecule', 'oxygen', 3/2)), 5.39),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'water', 3)), (('molecule', 'ethane', 1), ('molecule', 'oxygen', 5 / 2)), 8.94),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'water', 4)), (('molecule', 'propane', 1), ('molecule', 'oxygen', 7 / 2)), 12.39),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'water', 5)), (('molecule', 'butane', 1), ('molecule', 'oxygen', 9 / 2)), 15.8),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'water', 6)), (('molecule', 'pentane', 1), ('molecule', 'oxygen', 11 / 2)), 19.24),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', 'hexane', 1), ('molecule', 'oxygen', 13 / 2)), 22.69),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', 'heptane', 1), ('molecule', 'oxygen', 15 / 2)), 26.13),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'water', 9)), (('molecule', 'octane', 1), ('molecule', 'oxygen', 17 / 2)), 29.56),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'water', 11)), (('molecule', 'decane', 1), ('molecule', 'oxygen', 21 / 2)), 36.44),
)

iso_alkane_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'water', 5)), (('molecule', 'isobutane', 1), ('molecule', 'oxygen', 9 / 2)), 15.72),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'water', 6)), (('molecule', 'isopentane', 1), ('molecule', 'oxygen', 11 / 2)), 19.17),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', 'isohexane', 1), ('molecule', 'oxygen', 13 / 2)), 22.61),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', '3-methylpentane', 1), ('molecule', 'oxygen', 13 / 2)), 22.64),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', 'isoheptane', 1), ('molecule', 'oxygen', 15 / 2)), 26.05),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', '3-methylhexane', 1), ('molecule', 'oxygen', 15 / 2)), 26.09),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'water', 9)), (('molecule', '2-methylheptane', 1), ('molecule', 'oxygen', 17 / 2)), 29.49),
)

neo_alakane_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'water', 6)), (('molecule', 'neopentane', 1), ('molecule', 'oxygen', 11 / 2)), 19.02),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', 'neohexane', 1), ('molecule', 'oxygen', 13 / 2)), 22.49),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', 'neoheptane', 1), ('molecule', 'oxygen', 15 / 2)), 25.94),
)

amine_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'ammonia', 1), ('molecule', 'oxygen', 5 / 4)), 2.34),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'methylamine', 1), ('molecule', 'oxygen', 9 / 4)), 6.23),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'ethylamine', 1), ('molecule', 'oxygen', 13 / 4)), 9.62),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'propylamine', 1), ('molecule', 'oxygen', 17 / 4)), 13.04),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'isopropylamine', 1), ('molecule', 'oxygen', 17 / 4)), 12.90),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '1-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.47),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '2-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.34),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'isobutylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.4),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'tert-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.17),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'dimethylamine', 1), ('molecule', 'oxygen', 13 / 4)), 9.92),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'diethylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.67),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'methyltertbutylamine', 1), ('molecule', 'oxygen', 25 / 4)), 19.95),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'dipropylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.52),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'diisopropylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.23),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 17 / 2)), (('molecule', 'isopropyl-tertbutylamine', 1), ('molecule', 'oxygen', 33 / 4)), 26.67),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'dibutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.41),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'butylisobutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.26),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'diisobutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.17),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'ditertbutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.25),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'triethylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.76),
    reaction((('molecule', 'carbon-monoxide', 9), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 21 / 2)), (('molecule', 'tripropylamine', 1), ('molecule', 'oxygen', 41 / 4)), 34.01),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'trimethylamine', 1), ('molecule', 'oxygen', 17 / 4)), 13.52),
)

nitro_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'nitromethane', 1), ('molecule', 'oxygen', 3 / 4)), 3.12),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'nitroethane', 1), ('molecule', 'oxygen', 7 / 4)), 6.54),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '1-nitropropane', 1), ('molecule', 'oxygen', 11 / 4)), 9.97),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '2-nitropropane', 1), ('molecule', 'oxygen', 11 / 4)), 9.82),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '1-nitrobutane', 1), ('molecule', 'oxygen', 15 / 4)), 13.42),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-nitrobutane', 1), ('molecule', 'oxygen', 15 / 4)), 13.22),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '1-nitropentane', 1), ('molecule', 'oxygen', 19 / 4)), 16.86),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-methyl-2-nitropropane', 1), ('molecule', 'oxygen', 15 / 4)), 13.08),
)

nitrate_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitric-oxide', 1), ('molecule', 'water',  2)), (('molecule', 'nitric-acid', 1), ('molecule', 'hydrogen', 3 / 2)), 2.68),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water',  3 / 2)), (('molecule', 'methylnitrate', 1), ('molecule', 'oxygen', 1 / 4)), 2.69),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water',  5 / 2)), (('molecule', 'ethylnitrate', 1), ('molecule', 'oxygen', 5 / 4)), 6.01),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water',  7 / 2)), (('molecule', 'propylnitrate', 1), ('molecule', 'oxygen', 9 / 4)), 9.46),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water',  7 / 2)), (('molecule', 'isopropylnitrate', 1), ('molecule', 'oxygen', 9 / 4)), 9.28),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 3), ('molecule', 'water',  3)), (('molecule', '1.2.3-propanetriol-trinitrate', 1), ('molecule', 'hydrogen', 1 / 2)), 5.22),
)

nitrite_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitric-oxide', 1), ('molecule', 'water', 1)), (('molecule', 'nitrous-acid-trans', 1), ('molecule', 'hydrogen', 1 / 2)), 0.74),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'methyl-nitrite', 1), ('molecule', 'oxygen', 3 / 4)), 3.27),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'ethyl-nitrite', 1), ('molecule', 'oxygen', 7 / 4)), 6.49),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'propyl-nitrite', 1), ('molecule', 'oxygen', 11 / 4)), 10.03),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'isopropyl-nitrite', 1), ('molecule', 'oxygen', 11 / 4)), 9.88),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'butyl-nitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.4),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'sec-butyl-nitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.33),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'tert-butylnitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.13),
)

hydroxylamine_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'hydroxylamine', 1), ('molecule', 'oxygen', 3 / 4)), 3.31),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n.n-diethylhydroxylamine', 1), ('molecule', 'oxygen', 19 / 4)), 16.16),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'o-methylhydroxylamine', 1), ('molecule', 'oxygen', 7 / 4)), 6.21),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'n-methylhydroxylamine', 1), ('molecule', 'oxygen', 7 / 4)), 5.94),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'n.o-dimethylhydroxylamine', 1), ('molecule', 'oxygen', 11 / 4)), 9.73),
)

aromatic_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 3)), (('molecule', 'benzene', 1), ('molecule', 'oxygen', 9 / 2)), 15.25),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 4)), (('molecule', 'toluene', 1), ('molecule', 'oxygen', 11 / 2)), 18.57),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'water', 5)), (('molecule', 'ethylbenzene', 1), ('molecule', 'oxygen', 13 / 2)), 22.01),
    reaction((('molecule', 'carbon-monoxide', 9), ('molecule', 'water', 6)), (('molecule', 'propylbenzene', 1), ('molecule', 'oxygen', 15 / 2)), 25.43),
    reaction((('molecule', 'carbon-monoxide', 9), ('molecule', 'water', 6)), (('molecule', 'isopropylbenzene', 1), ('molecule', 'oxygen', 15 / 2)), 25.39),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'water', 7)), (('molecule', 'butylbenzene', 1), ('molecule', 'oxygen', 17 / 2)), 28.88),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'water', 7)), (('molecule', 'sec-butylbenzene', 1), ('molecule', 'oxygen', 17 / 2)), 28.81),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'water', 7)), (('molecule', 'tert-butylbenzene', 1), ('molecule', 'oxygen', 17 / 2)), 28.76),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'water', 7)), (('molecule', 'isobutylbenzene', 1), ('molecule', 'oxygen', 17 / 2)), 28.77),
)

aniline_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'aniline', 1), ('molecule', 'oxygen', 21 / 4)), 15.61),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-methylaniline', 1), ('molecule', 'oxygen', 25 / 4)), 18.94),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'n-methylaniline', 1), ('molecule', 'oxygen', 25 / 4)), 19.22),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n-ethylaniline', 1), ('molecule', 'oxygen', 29 / 4)), 22.59),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n.n-dimethylaniline', 1), ('molecule', 'oxygen', 29 / 4)), 23.05),
)

hydrazine_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'nitric-oxide', 2), ('molecule', 'water', 2)), (('molecule', 'hydrazine', 1), ('molecule', 'oxygen', 2)), 4.11),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 2), ('molecule', 'water', 4)), (('molecule', 'phenylhydrazine', 1), ('molecule', 'oxygen', 6)), 17.11),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 2), ('molecule', 'water', 4)), (('molecule', '1.2-dimethylhydrazine', 1), ('molecule', 'oxygen', 4)), 11.38),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 2), ('molecule', 'water', 4)), (('molecule', '1.1-dimethylhydrazine', 1), ('molecule', 'oxygen', 4)), 11.3),
)

amide_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'formamide', 1), ('molecule', 'oxygen', 5 / 4)), 1.95),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'acetamide', 1), ('molecule', 'oxygen', 9 / 4)), 5.14),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'propanamide', 1), ('molecule', 'oxygen', 13 / 4)), 8.58),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'butanamide', 1), ('molecule', 'oxygen', 17 / 4)), 11.99),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'pentanamide', 1), ('molecule', 'oxygen', 21 / 4)), 15.56),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'hexanamide', 1), ('molecule', 'oxygen', 25 / 4)), 18.86),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 17 / 2)), (('molecule', 'octanamide', 1), ('molecule', 'oxygen', 33 / 4)), 25.76),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'isobutiramide', 1), ('molecule', 'oxygen', 17 / 4)), 11.99),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '2.2-dimethylpropanamide', 1), ('molecule', 'oxygen', 21 / 4)), 15.32),
)

nitrile_gaseous: Tuple[reaction, ...] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 1 / 2)), (('molecule', 'hydrogen-cyanide', 1), ('molecule', 'oxygen', 5 / 4)), 2.85),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'acetonitrile', 1), ('molecule', 'oxygen', 9 / 4)), 5.87),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'propanenitrile', 1), ('molecule', 'oxygen', 13 / 4)), 9.29),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'butanenitrile', 1), ('molecule', 'oxygen', 17 / 4)), 12.76),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'pentanenitrile', 1), ('molecule', 'oxygen', 21 / 4)), 16.17),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'heptanenitrile', 1), ('molecule', 'oxygen', 29 / 4)), 23.04),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'octanenitrile', 1), ('molecule', 'oxygen', 33 / 4)), 26.49),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'decanenitrile', 1), ('molecule', 'oxygen', 41 / 4)), 33.37),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitric-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '2-methylpropanenitrile', 1), ('molecule', 'oxygen', 17 / 4)), 12.65),
)

all_gaseous_reactions = simple_alkane_gaseous + iso_alkane_gaseous + neo_alakane_gaseous + amine_gaseous + nitro_gaseous + nitrate_gaseous + nitrite_gaseous + hydroxylamine_gaseous + aromatic_gaseous + aniline_gaseous + hydrazine_gaseous + amide_gaseous + nitrile_gaseous
all_gaseous_reactions_named = {
    'Simple alkanes': simple_alkane_gaseous,
    'Branched alkanes (iso)': iso_alkane_gaseous,
    'Branched alkanes (neo)': neo_alakane_gaseous,
    'Amines': amine_gaseous,
    'Nitros': nitro_gaseous,
    'Nitrates': nitrate_gaseous,
    'Nitrites': nitrite_gaseous,
    'Hydroxylamines': hydroxylamine_gaseous,
    'Aromatics': aromatic_gaseous,
    'Anilines': aniline_gaseous,
    'Hydrazines': hydrazine_gaseous,
    'Amides': amide_gaseous,
    'Nitriles': nitrile_gaseous,
}