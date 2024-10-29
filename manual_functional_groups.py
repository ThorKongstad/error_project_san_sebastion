

molecule_functional_dict = {
# alkanes
    'methane': {'methyl_carbons': 1},
    'ethane': {'methyl_carbons': 2},
    'propane': {'methyl_carbons': 3},
    'butane': {'methyl_carbons': 4},
    'pentane': {'methyl_carbons': 5},
    'hexane': {'methyl_carbons': 6},
    'heptane': {'methyl_carbons': 7},
    'octane': {'methyl_carbons': 8},
    'decane': {'methyl_carbons': 9},

# iso-carbons
    'isobutane': {'iso_carbons': 1, 'methyl_carbons': 3},
    'isopentane': {'iso_carbons': 1, 'methyl_carbons': 4},
    'isohexane': {'iso_carbons': 1, 'methyl_carbons': 5},
    '3-methylpentane': {'iso_carbons': 1, 'methyl_carbons': 5},
    'isoheptane': {'iso_carbons': 1, 'methyl_carbons': 6},
    '3-methylhexane': {'iso_carbons': 1, 'methyl_carbons': 6},
    '2-methylheptane': {'iso_carbons': 1, 'methyl_carbons': 7},

# neo-carbons
    'neopentane': {'neo_carbons': 1, 'methyl_carbons': 4},
    'neohexane': {'neo_carbons': 1, 'methyl_carbons': 5},
    'neoheptane': {'neo_carbons': 1, 'methyl_carbons': 6},

# amines
    'ammonia': {'amines': 1},
    'methylamine': {'amines': 1, 'methyl_carbons': 1},
    'ethylamine': {'amines': 1, 'methyl_carbons': 2},
    'propylamine': {'amines': 1, 'methyl_carbons': 3},
    'isopropylamine': {'amines': 1, 'methyl_carbons': 2,'iso_carbons': 1},
    '1-butylamine': {'amines': 1, 'methyl_carbons': 4},
    '2-butylamine': {'amines': 1, 'methyl_carbons': 4},
    'isobutylamine': {'amines': 1, 'methyl_carbons': 3, 'iso_carbons': 1},
    'tert-butylamine': {'amines': 1, 'methyl_carbons': 3, 'iso_carbons': 1},
    'dimethylamine': {'amines': 1, 'methyl_carbons': 2},
    'diethylamine': {'amines': 1, 'methyl_carbons': 4},
    'methyltertbutylamine': {'amines': 1, 'methyl_carbons': 4, 'iso_carbons': 1},
    'dipropylamine': {'amines': 1, 'methyl_carbons': 6},
    'diisopropylamine': {'amines': 1, 'methyl_carbons': 6},
    'isopropyl-tertbutylamine': {'amines': 1, 'methyl_carbons': 6, 'iso_carbons': 1},
    'dibutylamine': {'amines': 1, 'methyl_carbons': 8},
    'butylisobutylamine': {'amines': 1, 'methyl_carbons': 7, 'iso_carbons': 1},
    'diisobutylamine': {'amines': 1, 'methyl_carbons': 6, 'iso_carbons': 2},
    'ditertbutylamine': {'amines': 1, 'methyl_carbons': 6, 'iso_carbons': 2},
    'trimethylamine': {'amines': 1, 'methyl_carbons': 1},
    'triethylamine': {'amines': 1, 'methyl_carbons': 6},
    'tripropylamine': {'amines': 1, 'methyl_carbons': 9},

# nitro
    'nitromethane': {'nitro': 1, 'methyl_carbons': 1},
    'nitroethane': {'nitro': 1, 'methyl_carbons': 2},
    '1-nitropropane': {'nitro': 1, 'methyl_carbons': 3},
    '2-nitropropane': {'nitro': 1, 'methyl_carbons': 3},
    '1-nitrobutane': {'nitro': 1, 'methyl_carbons': 4},
    '2-nitrobutane': {'nitro': 1, 'methyl_carbons': 4},
    '1-nitropentane': {'nitro': 1, 'methyl_carbons': 5},
    '2-methyl-2-nitropropane': {'nitro': 1, 'methyl_carbons': 3, 'iso_carbons': 1},

# nitrates
    'nitric-acid': {'nitrate': 1},
    'methylnitrate': {'nitrate': 1, 'methyl_carbons': 1},
    'ethylnitrate': {'nitrate': 1, 'methyl_carbons': 2},
    'propylnitrate': {'nitrate': 1, 'methyl_carbons': 3},
    'isopropylnitrate': {'nitrate': 1, 'methyl_carbons': 3},
    '1.2.3-propanetriol-trinitrate': {'nitrate': 3, 'methyl_carbons': 3},

#nitrites
    'nitrous-acid': {'nitrite': 1},
    'methyl-nitrite': {'nitrite': 1, 'methyl_carbons': 1},
    'ethyl-nitrite': {'nitrite': 1, 'methyl_carbons': 2},
    'propyl-nitrite': {'nitrite': 1, 'methyl_carbons': 3},
    'isopropyl-nitrite': {'nitrite': 1, 'methyl_carbons': 3},
    'butyl-nitrite': {'nitrite': 1, 'methyl_carbons': 4},
    'sec-butyl-nitrite': {'nitrite': 1, 'methyl_carbons': 4},
    'tert-butylnitrite': {'nitrite': 1, 'methyl_carbons': 3, 'iso_carbons': 1},

# hydroxylamine
    'hydroxylamine': {'hydroxylamine': 1},
    'n.n-diethylhydroxylamine': {'hydroxylamine': 1, 'methyl_carbons': 4},
    'o-methylhydroxylamine': {'hydroxylamine': 1, 'methyl_carbons': 1},
    'n-methylhydroxylamine': {'hydroxylamine': 1, 'methyl_carbons': 1},
    'n.o-dimethylhydroxylamine': {'hydroxylamine': 1, 'methyl_carbons': 2},

# aromatics
    'benzene': {'phenyl': 1},
    'toluene': {'phenyl': 1, 'methyl_carbons': 1},
    'ethylbenzene': {'phenyl': 1, 'methyl_carbons': 2},
    'propylbenzene': {'phenyl': 1, 'methyl_carbons': 3},
    'isopropylbenzene': {'phenyl': 1, 'methyl_carbons': 2, 'iso_carbons': 1},
    'butylbenzene': {'phenyl': 1, 'methyl_carbons': 4},
    'sec-butylbenzene': {'phenyl': 1, 'methyl_carbons': 3, 'iso_carbons': 1},
    'tert-butylbenzene': {'phenyl': 1, 'methyl_carbons': 3, 'neo_carbons': 1},
    'isobutylbenzene': {'phenyl': 1, 'methyl_carbons': 3, 'iso_carbons': 1},

# aniline
    'aniline': {'aniline': 1},
    '2-methylaniline': {'aniline': 1, 'methyl_carbons': 1},
    'n-methylaniline': {'aniline': 1, 'methyl_carbons': 1},
    'n-ethylaniline': {'aniline': 1, 'methyl_carbons': 2},
    'n.n-dimethylaniline': {'aniline': 1, 'methyl_carbons': 2},

# hydrazine
    'phenylhydrazine': {'phenyl': 1, 'hydrazine': 1},
    '1.2-dimethylhydrazine': {'hydrazine': 1, 'methyl_carbons': 2},
    '1.1-dimethylhydrazine': {'hydrazine': 1, 'methyl_carbons': 2},
    'hydrazine': {'hydrazine': 1},

# amide
    'formamide': {'amide': 1, 'methyl_carbons': 1},
    'acetamide': {'amide': 1, 'methyl_carbons': 2},
    'propanamide': {'amide': 1, 'methyl_carbons': 3},
    'butanamide': {'amide': 1, 'methyl_carbons': 4},
    'pentanamide': {'amide': 1, 'methyl_carbons': 5},
    'hezanamide': {'amide': 1, 'methyl_carbons': 6},
    'octanamide': {'amide': 1, 'methyl_carbons': 8},
    'isobutiramide': {'amide': 1, 'methyl_carbons': 3, 'iso_carbons': 1},
    '2.2-dimethylpropanamide': {'amide': 1, 'methyl_carbons': 3, 'neo_carbons': 1},

# nitrile
    'hydrogencyanide': {'nitrile': 1},
    'acetonitrile': {'nitrile': 1, 'methyl_carbons': 1},
    'propanenitrile': {'nitrile': 1, 'methyl_carbons': 2},
    'butanenitrile': {'nitrile': 1, 'methyl_carbons': 3},
    'pentanenitrile': {'nitrile': 1, 'methyl_carbons': 4},
    'heptanenitrile': {'nitrile': 1, 'methyl_carbons': 6},
    'octanenitrile': {'nitrile': 1, 'methyl_carbons': 7},
    'decanenitrile': {'nitrile': 1, 'methyl_carbons': 9},
    '2-methylpropanenitrile': {'nitrile': 1, 'methyl_carbons': 2, 'iso_carbons': 1},
}