

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
    'decane': {'methyl_carbons': 10}, # was 9

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
    'ammonia': {'amines': 1}, # Is an amine but do not correct it. It was assumed to be well described and therefore is the basis for N2 error calculation
    'methylamine': {'amines': 1, 'methyl_carbons': 1},
    'ethylamine': {'amines': 1, 'methyl_carbons': 2},
    'propylamine': {'amines': 1, 'methyl_carbons': 3},
    'isopropylamine': {'amines': 1, 'methyl_carbons': 2, 'iso_carbons': 1},
    '1-butylamine': {'amines': 1, 'methyl_carbons': 4},
    '2-butylamine': {'amines': 1, 'methyl_carbons': 3, 'iso_carbons': 1}, # was 4 methyl_carbon
    'isobutylamine': {'amines': 1, 'methyl_carbons': 3, 'iso_carbons': 1}, 
    'tert-butylamine': {'amines': 1, 'methyl_carbons': 3, 'neo_carbons': 1}, # was iso_carbon instaed of neo_carbons
    'dimethylamine': {'amines': 1, 'methyl_carbons': 2},
    'diethylamine': {'amines': 1, 'methyl_carbons': 4},
    'methyltertbutylamine': {'amines': 1, 'methyl_carbons': 4, 'neo_carbons': 1}, # was iso_carbon instaed of neo_carbons
    'dipropylamine': {'amines': 1, 'methyl_carbons': 6}, 
    'diisopropylamine': {'amines': 1, 'methyl_carbons': 4, 'iso_carbons': 2}, # was 6 methyl_carbons
    'isopropyl-tertbutylamine': {'amines': 1, 'methyl_carbons': 5, 'iso_carbons': 1, 'neo_carbons': 1},  # was 6 methyl_carbons and no neo_carbon
    'dibutylamine': {'amines': 1, 'methyl_carbons': 8},
    'butylisobutylamine': {'amines': 1, 'methyl_carbons': 7, 'iso_carbons': 1},
    'diisobutylamine': {'amines': 1, 'methyl_carbons': 6, 'iso_carbons': 2},
    'ditertbutylamine': {'amines': 1, 'methyl_carbons': 6, 'neo_carbons': 2}, # was 2 iso_carbon
    'trimethylamine': {'amines': 1, 'methyl_carbons': 3}, # was 1 methyl_carbons 
    'triethylamine': {'amines': 1, 'methyl_carbons': 6},
    'tripropylamine': {'amines': 1, 'methyl_carbons': 9},

# nitro
    'nitromethane': {'nitro': 1, 'methyl_carbons': 1},
    'nitroethane': {'nitro': 1, 'methyl_carbons': 2},
    '1-nitropropane': {'nitro': 1, 'methyl_carbons': 3},
    '2-nitropropane': {'nitro': 1, 'methyl_carbons': 2, 'iso_carbons': 1}, # was 3 methyl_carbons
    '1-nitrobutane': {'nitro': 1, 'methyl_carbons': 4},
    '2-nitrobutane': {'nitro': 1, 'methyl_carbons': 3, 'iso_carbons': 1}, # was 4 methyl_carbons
    '1-nitropentane': {'nitro': 1, 'methyl_carbons': 5},
    '2-methyl-2-nitropropane': {'nitro': 1, 'methyl_carbons': 3, 'neo_carbons': 1}, # was 1 iso_carbon 

# nitrates
    'nitric-acid': {'nitrate': 1},
    'methylnitrate': {'nitrate': 1, 'methyl_carbons': 1},
    'ethylnitrate': {'nitrate': 1, 'methyl_carbons': 2},
    'propylnitrate': {'nitrate': 1, 'methyl_carbons': 3},
    'isopropylnitrate': {'nitrate': 1, 'methyl_carbons': 2, 'iso_carbons': 1}, # was 3 methyl_carbons
    '1.2.3-propanetriol-trinitrate': {'nitrate': 3, 'methyl_carbons': 2, 'iso_carbons': 1}, # was 3 methyl_carbons

#nitrites
    'nitrous-acid': {'nitrite': 1},
    'methyl-nitrite': {'nitrite': 1, 'methyl_carbons': 1},
    'ethyl-nitrite': {'nitrite': 1, 'methyl_carbons': 2},
    'propyl-nitrite': {'nitrite': 1, 'methyl_carbons': 3},
    'isopropyl-nitrite': {'nitrite': 1, 'methyl_carbons': 2, 'iso_carbons': 1}, # was 3 methyl_carbons
    'butyl-nitrite': {'nitrite': 1, 'methyl_carbons': 4},
    'sec-butyl-nitrite': {'nitrite': 1, 'methyl_carbons': 3, 'iso_carbons': 1}, # was 4 methyl_carbons
    'tert-butylnitrite': {'nitrite': 1, 'methyl_carbons': 3, 'neo_carbons': 1}, # was 1 iso_carbons

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
    'formamide': {'amide': 1}, # was 1 methyl_carbons
    'acetamide': {'amide': 1, 'methyl_carbons': 1}, # was 2 methyl_carbons
    'propanamide': {'amide': 1, 'methyl_carbons': 2}, # was 3 methyl_carbons
    'butanamide': {'amide': 1, 'methyl_carbons': 3},  # was 4 methyl_carbons
    'pentanamide': {'amide': 1, 'methyl_carbons': 4}, # was 5 methyl_carbons
    'hexanamide': {'amide': 1, 'methyl_carbons': 5}, # was hezanadamide and 6 methyl_carbons 
    'octanamide': {'amide': 1, 'methyl_carbons': 7}, # was 8 methyl_carbons 
    'isobutiramide': {'amide': 1, 'methyl_carbons': 2, 'iso_carbons': 1}, # was 3 methyl_carbons
    '2.2-dimethylpropanamide': {'amide': 1, 'methyl_carbons': 3, 'neo_carbons': 1},

# nitrile
    'hydrogen-cyanide': {'nitrile': 1},
    'acetonitrile': {'nitrile': 1, 'methyl_carbons': 1},
    'propanenitrile': {'nitrile': 1, 'methyl_carbons': 2},
    'butanenitrile': {'nitrile': 1, 'methyl_carbons': 3},
    'pentanenitrile': {'nitrile': 1, 'methyl_carbons': 4},
    'heptanenitrile': {'nitrile': 1, 'methyl_carbons': 6},
    'octanenitrile': {'nitrile': 1, 'methyl_carbons': 7},
    'decanenitrile': {'nitrile': 1, 'methyl_carbons': 9},
    '2-methylpropanenitrile': {'nitrile': 1, 'methyl_carbons': 2, 'iso_carbons': 1},
}