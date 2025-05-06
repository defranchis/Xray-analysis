
sample_tags = {
    "1006_LR":  {"tag": "FZ -5V",                   "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1113_LR":  {"tag": "FZ -2V",                   "structures": ["MOShalf", "MOS2000", None]},
    "3001_UL":  {"tag": "EPI -5V",                  "structures": ["MOShalf", "MOS2000", "GCD"]},
    "3007_UL":  {"tag": "EPI -2V",                  "structures": [None, "MOS2000", "GCD"]},
    "3103_LR":  {"tag": "EPI -2V",                  "structures": ["MOShalf", None, None]},

    "1008_LR":  {"tag": "Type B -5V",               "structures": ["MOShalf", "MOS2000", None]},
    "1009_LR":  {"tag": "Type C -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1010_UL":  {"tag": "Type D -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},
    "1011_LR":  {"tag": "Type E -5V",               "structures": ["MOShalf", "MOS2000", "GCD"]},

    "3010_LR":  {"tag": "HGCAL 6inch -5V",          "structures": ["MOShalf", "MOS2000", "GCD"]},
    "23_SE_GCD": {"tag": "Tracker -2V",             "structures": [None, "MOS", "GCD"]},
    "24_E_MOS": {"tag": "Tracker -2V",              "structures": ["MOSc1", None, None]},
    "N0538_25_LR": {"tag": "New Type C -2V",        "structures": ["MOShalf", "MOS2000", "GCD"]}
}

    #"1003_LR":  {"tag": "FZ -2V (rep.)",            "structures": ["MOShalf", "MOS2000", "GCD"]},
    #"1012_UL":  {"tag": "FZ -2V (rep.)",            "structures": [None, None, "GCD"]},
    #"1105_LR":  {"tag": "FZ -5V (rep.)",            "structures": ["MOShalf", "MOS2000", None]},
    #"1109_LR":  {"tag": "Type B -5V (rep.)",        "structures": ["MOShalf", "MOS2000", None]},
    #"1112_LR":  {"tag": "Type E -5V (rep.)",        "structures": ["MOShalf", "MOS2000", None]},
    #"3003_UL":  {"tag": "EPI -2V (rep.)",           "structures": ["MOShalf", "MOS2000", "GCD"]},
    #"3009_LR":  {"tag": "HGCAL 6inch -5V (rep.)",   "structures": ["MOShalf", "MOS2000", "GCD"]},
    #"3101_LR":  {"tag": "EPI -5V (rep.)",           "structures": ["MOShalf", "MOS2000", "GCD"]},
    #"24_E_MOS": {"tag": "Tracker square",           "structures": ["MOSs2", None, None]},
    #"24_E_MOS": {"tag": "Tracker circular (2)",     "structures": [None, "MOSc2", None]},

tag_styles = {
    "FZ -5V":               {"color": "black",      "marker": "o"},
    "FZ -2V":               {"color": "magenta",    "marker": "D"},
    "EPI -5V":              {"color": "darkgreen",  "marker": "^"},
    "EPI -2V":              {"color": "orange",     "marker": "v"},
    "Type B -5V":           {"color": "blue",       "marker": "^"},
    "Type C -5V":           {"color": "red",        "marker": "v"},
    "Type D -5V":           {"color": "gold",       "marker": "s"},
    "Type E -5V":           {"color": "green",      "marker": ">"},
    "HGCAL 6inch -5V":      {"color": "blue",       "marker": "v"},
    "Tracker -2V":          {"color": "cyan",       "marker": "^"},
    "New Type C -2V":       {"color": "darkorange", "marker": "D"},
}

def getAllFloatingMOS():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][0] is not None]

def getAllBiasedMOS():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][1] is not None]

def getAllGCD():
    return [x for x in sample_tags.keys() if sample_tags[x]["structures"][2] is not None]