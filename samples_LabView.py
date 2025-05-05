
sample_tags = {
    "1006_LR":  {"tag": "Standard -5V",             "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1008_LR":  {"tag": "Type B -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1009_LR":  {"tag": "Type C -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1010_UL":  {"tag": "Type D -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1011_LR":  {"tag": "Type E -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1003_LR":  {"tag": "Standard -2V",             "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1113_LR":  {"tag": "Standard -2V",             "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3009_LR":  {"tag": "HGCAL 6'' -5V",            "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3001_UL":  {"tag": "EPI -5V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1112_LR":  {"tag": "Type E -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3003_UL":  {"tag": "EPI -2V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3103_LR":  {"tag": "EPI -2V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1109_LR":  {"tag": "Type B -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1105_LR":  {"tag": "Standard -5V",             "structures": ["MOShalf", "MOS2000", None]},
    "3101_LR":  {"tag": "EPI -5V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3010_LR":  {"tag": "HGCAL 6'' -5V",            "structures": ["MOShalf", "MOS2000", "GCD"]},
    "24_E_MOS": {"tag": "Tracker",                  "structures": ["MOSc1", "MOSc2", None]},
    #"24_E_MOS": {"tag": "Tracker square",           "structures": ["MOSs2", None, None]},
    #"24_E_MOS": {"tag": "Tracker circular (2)",     "structures": [None, "MOSc2", None]},
    "N0538_25_LR":     {"tag": "New Type C",        "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3007_UL":  {"tag": "EPI -2V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1012_UL":  {"tag": "FZ -2V (GDC repetition)",  "structures": [None, None, "GCD"]},
    "23_SE_GCD": {"tag": "tracker biased",          "structures": [None, "MOS", "GCD"]}
}

def getAllFloatingMOS():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][0] is not None]

def getAllBiasedMOS():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][1] is not None]

def getAllGCD():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][2] is not None]