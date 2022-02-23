## safe batch style
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# # leg= legend entry text
# # lc = linecolor
# # lw = linewidth
# # ls = linestyle
# # ms = markerstyle
sampleNames = {
    'N4791-1_LR'  : {'leg': 'A FZ',             'lw'  : 2, 'ls'  : 1, 'lc'  : ROOT.kBlack,    'ms' : ROOT.kFullCircle},
    'N4790-1_UL'  : {'leg': 'B FZ',             'lw'  : 2, 'ls'  : 1, 'lc'  : ROOT.kGreen+2,  'ms' : ROOT.kOpenCircle},
    'N4790-1_LR'  : {'leg': 'B FZ #2',          'lw'  : 2, 'ls'  : 3, 'lc'  : ROOT.kGreen+1,    'ms' : ROOT.kOpenTriangleDown},
    'N4791-6_UL'  : {'leg': 'C FZ',             'lw'  : 2, 'ls'  : 1, 'lc'  : ROOT.kAzure+1,      'ms' : ROOT.kFullTriangleDown},
    'N4791-6_LR'  : {'leg': 'C FZ #2',          'lw'  : 2, 'ls'  : 1, 'lc'  : ROOT.kRed+1,     'ms' : ROOT.kOpenTriangleDown}, # or kBlue when comparing multiple type C together
    'N4789-10_UL' : {'leg': 'C EPI',   'lw'  : 2, 'ls'  : 7, 'lc'  : ROOT.kRed+2,    'ms' : ROOT.kFullSquare},
    #'N4789-10_UL' : {'leg': 'C EPI 14 kGy/h (reference)',   'lw'  : 2, 'ls'  : 7, 'lc'  : ROOT.kRed+2,    'ms' : ROOT.kFullSquare},
    'N4789-10_LR' : {'leg': 'C EPI #2',         'lw'  : 2, 'ls'  : 3, 'lc'  : ROOT.kRed,    'ms' : ROOT.kOpenSquare},
    'N4789-12_UL' : {'leg': 'C EPI 39 kGy/h',   'lw'  : 2, 'ls'  : 7, 'lc'  : ROOT.kOrange+2, 'ms' : ROOT.kOpenTriangleUp},
    'N4789-12_LR' : {'leg': 'C EPI  7 kGy/h',   'lw'  : 2, 'ls'  : 7, 'lc'  : ROOT.kGray+3,   'ms' : ROOT.kOpenSquareDiagonal},
    'N4790-13_UL' : {'leg': 'D FZ',             'lw'  : 2, 'ls'  : 1, 'lc'  : ROOT.kAzure,    'ms' : ROOT.kFullCircle},
    'N4790-13_LR' : {'leg': 'D FZ #2',          'lw'  : 2, 'ls'  : 3, 'lc'  : ROOT.kAzure+2,  'ms' : ROOT.kOpenCircle},
    'N4788-9_LR'  : {'leg': 'D EPI',            'lw'  : 2, 'ls'  : 7, 'lc'  : ROOT.kBlue,     'ms' : ROOT.kOpenSquare},
}

def getSampleNames():
    return [x for x in sampleNames.keys()]

def getSampleAttribute(name, key):
    return sampleNames[name][key];
