import sys
import pathlib
from typing import Tuple

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
from error_project_san_sebastion.reaction_functions import reaction

#reaction((('',),), (('',),), experimental value)
# https://doi.org/10.1002/cctc.202100125

simple_alkane_formation: Tuple[reaction] = (
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


iso_alkane_formation: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'water', 5)), (('molecule', 'isobutane', 1), ('molecule', 'oxygen', 9 / 2)), 15.72),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'water', 6)), (('molecule', 'isopentane', 1), ('molecule', 'oxygen', 11 / 2)), 19.17),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', 'isohexane', 1), ('molecule', 'oxygen', 13 / 2)), 22.61),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', '3-methylpentane', 1), ('molecule', 'oxygen', 13 / 2)), 22.64),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', 'isoheptane', 1), ('molecule', 'oxygen', 15 / 2)), 26.05),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', '3-methylhexane', 1), ('molecule', 'oxygen', 15 / 2)), 26.09),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'water', 9)), (('molecule', 'isoheptane', 1), ('molecule', 'oxygen', 17 / 2)), 26.09),
)

neo_alakane_formation: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'water', 6)), (('molecule', 'neopentane', 1), ('molecule', 'oxygen', 11 / 2)), 19.02),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'water', 7)), (('molecule', 'neohexane', 1), ('molecule', 'oxygen', 13 / 2)), 22.49),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'water', 8)), (('molecule', 'neoheptane', 1), ('molecule', 'oxygen', 15 / 2)), 25.94),

)

amine_formation: Tuple[reaction] = (
    reaction((('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'ammonia', 1), ('molecule', 'oxygen', 5 / 4)), 2.34),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'methylamine', 1), ('molecule', 'oxygen', 9 / 4)), 6.23),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'ethylamine', 1), ('molecule', 'oxygen', 13 / 4)), 9.62),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'propylamine', 1), ('molecule', 'oxygen', 17 / 4)), 13.04),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'isopropylamine', 1), ('molecule', 'oxygen', 17 / 4)), 12.90),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '1-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.47),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '2-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.34),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'isobutylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.4),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'tert-butylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.17),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'dimethylamine', 1), ('molecule', 'oxygen', 13 / 4)), 9.92),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'diethylamine', 1), ('molecule', 'oxygen', 21 / 4)), 16.67),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'methyl-tert-butylamine', 1), ('molecule', 'oxygen', 25 / 4)), 19.95),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'dipropylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.52),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'diisoprpylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.23),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 17 / 2)), (('molecule', 'isopropyl-terbutylamine', 1), ('molecule', 'oxygen', 33 / 4)), 26.67),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'dibutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.41),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'butylisobutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.26),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'diisobutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.17),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'diterbutylamine', 1), ('molecule', 'oxygen', 37 / 4)), 30.25),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'triethylamine', 1), ('molecule', 'oxygen', 29 / 4)), 23.76),
    reaction((('molecule', 'carbon-monoxide', 9), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 21 / 2)), (('molecule', 'tripropylamine', 1), ('molecule', 'oxygen', 41 / 4)), 34.01),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'trimethylamine', 1), ('molecule', 'oxygen', 17 / 4)), 13.52),
)

nitro_formation: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'nitromethane', 1), ('molecule', 'oxygen', 3 / 4)), 3.12),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'nitroethane', 1), ('molecule', 'oxygen', 7 / 4)), 6.54),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '1-nitropropane', 1), ('molecule', 'oxygen', 11 / 4)), 9.97),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '2-nitropropane', 1), ('molecule', 'oxygen', 11 / 4)), 9.82),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '1-nitrobutane', 1), ('molecule', 'oxygen', 15 / 4)), 13.42),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-nitrobutane', 1), ('molecule', 'oxygen', 15 / 4)), 13.22),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '1-nitrobutane', 1), ('molecule', 'oxygen', 19 / 4)), 16.86),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-methyl-2-nitropropane', 1), ('molecule', 'oxygen', 15 / 4)), 13.08),
)

nitrate_formation: Tuple[reaction] = (
    reaction((('molecule', 'nitrc-oxide', 1), ('molecule', 'water',  2)), (('molecule', 'nitric-acid', 1), ('molecule', 'hydrogen', 3 / 2)), 2.68),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water',  3 / 2)), (('molecule', 'methylnitrate', 1), ('molecule', 'oxygen', 1 / 4)), 2.69),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water',  5 / 2)), (('molecule', 'ethylnitrate', 1), ('molecule', 'oxygen', 5 / 4)), 6.01),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water',  7 / 2)), (('molecule', 'propylnitrate', 1), ('molecule', 'oxygen', 9 / 4)), 9.46),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water',  7 / 2)), (('molecule', 'isopropylnitrate', 1), ('molecule', 'oxygen', 9 / 4)), 9.28),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 3), ('molecule', 'water',  3)), (('molecule', '1.2.3-propanetriol-trinitrate', 1), ('molecule', 'hydrogen', 1 / 2)), 5.22),
)

nitrite_formation: Tuple[reaction] = (
    reaction((('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 1)), (('molecule', 'nitrous-acid', 1), ('molecule', 'hydrogen', 1 / 2)), 0.74),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'methylnitrite', 1), ('molecule', 'oxygen', 3 / 4)), 3.27),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'ethylnitrite', 1), ('molecule', 'oxygen', 7 / 4)), 6.49),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'propylnitrite', 1), ('molecule', 'oxygen', 11 / 4)), 10.03),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'isopropylnitrite', 1), ('molecule', 'oxygen', 11 / 4)), 9.88),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'butyl-nitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.4),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'sec-butyl-nitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.33),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'tert-butyl-nitrite', 1), ('molecule', 'oxygen', 15 / 4)), 13.13),
)

hydroxylamineamine_formation: Tuple[reaction] = (
    reaction((('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'hydroxylamine', 1), ('molecule', 'oxygen', 3 / 4)), 3.31),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n.n-diethylhydroxylamine', 1), ('molecule', 'oxygen', 19 / 4)), 16.16),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'o-methylhydroxylamine', 1), ('molecule', 'oxygen', 7 / 4)), 6.21),
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'n-methylhydroxylamine', 1), ('molecule', 'oxygen', 7 / 4)), 5.94),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'n.o-diethylhydroxylamine', 1), ('molecule', 'oxygen', 11 / 4)), 9.73),
)

aromatic_formation: Tuple[reaction] = (
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

aniline_formations: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'aniline', 1), ('molecule', 'oxygen', 21 / 4)), 15.61),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', '2-methylaniline', 1), ('molecule', 'oxygen', 25 / 4)), 18.94),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'n-methylaniline', 1), ('molecule', 'oxygen', 25 / 4)), 19.22),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n-ethylaniline', 1), ('molecule', 'oxygen', 29 / 4)), 22.59),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'n.n-dimethylaniline', 1), ('molecule', 'oxygen', 29 / 4)), 23.05),
)

hydrazine_formation: Tuple[reaction] = (
    reaction((('molecule', 'nitrc-oxide', 2), ('molecule', 'water', 2)), (('molecule', 'hydrazine', 1), ('molecule', 'oxygen', 2)), 4.11),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 2), ('molecule', 'water', 4)), (('molecule', 'phenyl-hydrazine', 1), ('molecule', 'oxygen', 6)), 17.11),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 2), ('molecule', 'water', 4)), (('molecule', '1.2-dimethylhydrazine', 1), ('molecule', 'oxygen', 4)), 11.38),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 2), ('molecule', 'water', 4)), (('molecule', '1.1-dimethylhydrazine', 1), ('molecule', 'oxygen', 4)), 11.3),
)

amide_formation: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'formamide', 1), ('molecule', 'oxygen', 5 / 4)), 1.95),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'acetamide', 1), ('molecule', 'oxygen', 9 / 4)), 5.14),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'propanamide', 1), ('molecule', 'oxygen', 13 / 4)), 8.58),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'butanamide', 1), ('molecule', 'oxygen', 17 / 4)), 11.99),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', 'pentanamide', 1), ('molecule', 'oxygen', 21 / 4)), 15.56),
    reaction((('molecule', 'carbon-monoxide', 6), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'hexanamide', 1), ('molecule', 'oxygen', 25 / 4)), 18.86),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 17 / 2)), (('molecule', 'octanamide', 1), ('molecule', 'oxygen', 33 / 4)), 25.76),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'isobutiramide', 1), ('molecule', 'oxygen', 17 / 4)), 11.99),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 11 / 2)), (('molecule', '2.2-dimethylpropanamide', 1), ('molecule', 'oxygen', 21 / 4)), 15.32),
)

nitrile_formation: Tuple[reaction] = (
    reaction((('molecule', 'carbon-monoxide', 1), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 1 / 2)), (('molecule', 'hydrogen-cyanide', 1), ('molecule', 'oxygen', 5 / 4)), 2.85),
    reaction((('molecule', 'carbon-monoxide', 2), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 3 / 2)), (('molecule', 'acetonitrile', 1), ('molecule', 'oxygen', 9 / 4)), 5.87),
    reaction((('molecule', 'carbon-monoxide', 3), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 5 / 2)), (('molecule', 'propanenitrile', 1), ('molecule', 'oxygen', 13 / 4)), 9.29),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', 'butanenitrile', 1), ('molecule', 'oxygen', 17 / 4)), 12.76),
    reaction((('molecule', 'carbon-monoxide', 5), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 9 / 2)), (('molecule', 'pentanenitrile', 1), ('molecule', 'oxygen', 21 / 4)), 16.17),
    reaction((('molecule', 'carbon-monoxide', 7), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 13 / 2)), (('molecule', 'heptanenitrile', 1), ('molecule', 'oxygen', 29 / 4)), 23.04),
    reaction((('molecule', 'carbon-monoxide', 8), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 15 / 2)), (('molecule', 'octanenitrile', 1), ('molecule', 'oxygen', 33 / 4)), 26.49),
    reaction((('molecule', 'carbon-monoxide', 10), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 19 / 2)), (('molecule', 'decanenitrile', 1), ('molecule', 'oxygen', 41 / 4)), 33.37),
    reaction((('molecule', 'carbon-monoxide', 4), ('molecule', 'nitrc-oxide', 1), ('molecule', 'water', 7 / 2)), (('molecule', '2-ethylpropanenitrile', 1), ('molecule', 'oxygen', 17 / 4)), 12.65),
)

all_reactions = simple_alkane_formation + iso_alkane_formation + neo_alakane_formation + amine_formation + nitro_formation + nitrate_formation + nitrite_formation + hydroxylamineamine_formation + aromatic_formation + aniline_formations + hydrazine_formation + amide_formation + nitrile_formation