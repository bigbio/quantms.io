"""
Constants used in quantmsio for transformation of data from different sources into the quantms.io format.
"""

TMT_CHANNELS = {
    "TMT10": [
        "TMT126",
        "TMT127C",
        "TMT127N",
        "TMT128C",
        "TMT128N",
        "TMT129C",
        "TMT129N",
        "TMT130C",
        "TMT130N",
        "TMT131",
    ],
    "TMT11": [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
    ],
    "TMT16": [
        "TMT126",
        "TMT127N",
        "TMT127C",
        "TMT128N",
        "TMT128C",
        "TMT129N",
        "TMT129C",
        "TMT130N",
        "TMT130C",
        "TMT131N",
        "TMT131C",
        "TMT132N",
        "TMT132C",
        "TMT133N",
        "TMT133C",
        "TMT134N",
    ],
    "TMT6": ["TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131"],
}

ITRAQ_CHANNEL = {
    "ITRAQ4": ["ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117"],
    "ITRAQ8": [
        "ITRAQ113",
        "ITRAQ114",
        "ITRAQ115",
        "ITRAQ116",
        "ITRAQ117",
        "ITRAQ118",
        "ITRAQ119",
        "ITRAQ121",
    ],
    # NO EXAMPLES.
}

SINGLE_PROTEIN = "single_protein"
GROUP_PROTEIN = "indistinguishable_protein_group"
PROTEIN_DETAILS = "protein_details"
