import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.CMS)

from copy import *
from utility import *

from glob import glob
import os
import constants as cnst

from functionsAnnealing import *

import argparse
import numpy as np


# Default definitions
indir = '/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_MOS_GCD/LabView_Setup/'
outdir = 'allplots'
outfiles = 'Vfb_files'
doses = [0,1,2,5,10,15,20,30,40,70,100,181,200,394,436,509,749,762,1030]
samples = ['1006_LR','1008_LR','1009_LR','1010_UL','1011_LR','1003_LR','1113_LR','3009_LR',
           '3001_UL','1112_LR','3003_UL','3103_LR','1109_LR','1105_LR','3101_LR',
           '3010_LR','24_E_MOS','23_SE_GCD','N0538_25_LR','3007_UL','1012_UL']
GCD_exclude = ['1008_LR','1113_LR','1105_LR','1112_LR']
MOS_exclude = ['1012_UL']

def checkGoodList(args, ls):
    for i,l1 in enumerate(ls):
        for j,l2 in enumerate(ls):
            if i==j: continue
            if l1 in l2: 
                if l1+'_1kHz' in l2: continue
                if 'GCD' in l1 and l1.split('_')[1]+'_'+l1.split('_')[2] in args.GCD_exclude:
                    continue
                print('\nWARNING: potential duplicate:', l1, l2)
                print('please check\n')
    return

def readGoodList(args):
    global good_GCD
    global good_MOS
    f_GCD = open(f'{args.configdir}/good_GCD.txt','r') 
    f_MOS = open(f'{args.configdir}/good_MOS.txt','r') 
    good_GCD = f_GCD.read().splitlines()
    good_MOS = f_MOS.read().splitlines()
    return

def readMOSexcludeDose(args):
    global MOS_exclude_dose
    f = open(f'{args.configdir}/MOS_exclude_dose.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    MOS_exclude_dose = [l.replace(' ','').split('&') for l in lines]
    for l in MOS_exclude_dose:
        if len(l) != 3:
            print(f"ERROR in file MOS_exclude_dose.txt: wrong number of arguments in entry below\n {l}")
            sys.exit()
        if l[1] != 'MOShalf' and l[1] != 'MOS2000':
            print(f"ERROR in file MOS_exclude_dose.txt: structure in entry below not recognised\n {l}")
            sys.exit()
    return

def isMOSexcluded(sample,structure,dose):
    for l in MOS_exclude_dose:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            return True
    return False


def readGCDexcludeDose(args):
    global GCD_exclude_dose
    f = open(f'{args.configdir}/GCD_exclude_dose.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    GCD_exclude_dose = [l.replace(' ','').split('&') for l in lines]
    for l in GCD_exclude_dose:
        if len(l) != 2:
            print(f"ERROR in file GCD_exclude_dose.txt: wrong number of arguments in entry below\n {f}")
            sys.exit()
    return

def isGCDexcluded(sample,dose):
    for l in GCD_exclude_dose:
        if l[0] == sample:
            if '+' in l[1]:
                d = float(l[1].replace('+',''))
                if dose >= d:
                    return True
            else:
                if dose == float(l[1]):
                    return True
    return False

def readRemovePoints(args):
    global remove_points
    f = open(f'{args.configdir}/remove_last_point.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    remove_points = [l.replace(' ','').split('&') for l in lines]
    for l in remove_points:
        if len(l) != 3 and len(l) != 4:
            print(f"ERROR in file remove_last_point.txt: wrong number of arguments in entry below\n {l}")
            sys.exit()
    return

def pointsToRemove(sample,structure,dose):
    for l in remove_points:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) == 4:
                return int(l[3])
            else:
                return 1
    return 0

def readGCDcuts(args):
    global GCD_cuts
    f = open(f'{args.configdir}/cut_GCD_curve.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    GCD_cuts = [l.replace(' ','').split('&') for l in lines]
    for l in GCD_cuts:
        if len(l) != 4 and len(l) != 6:
            print("ERROR in file cut_GCD_curve.txt: wrong number of arguments in entry below\n {l}")
            sys.exit()
    return

def getGCDcuts(sample,dose):
    low = high = gmin = gmax = -999
    for l in GCD_cuts:
        if sample != l[0]: continue
        if '+' in l[1]:
            d = float(l[1].replace('+',''))
            if dose < d: continue
        elif '-' in l[1]:
            d = float(l[1].replace('-',''))
            if dose > d: continue
        else:
            d = float(l[1])
            if dose != d: continue
        low = float(l[2])
        high = float(l[3])
        if len(l) == 6:
            gmin = float(l[4])
            gmax = float(l[5])

    return low, high, gmin, gmax

def readGCDranges(args):
    global GCD_ranges
    f = open(f'{args.configdir}/GCD_range.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    GCD_ranges = [l.replace(' ','').split('&') for l in lines]
    for l in GCD_ranges:
        if len(l) != 4:
            print("ERROR in file GCD_range.txt: wrong number of arguments in entry below\n {l}")
            sys.exit()
    return

def getGCDrange(sample,dose,low,high):
    for l in GCD_ranges:
        if sample != l[0]: continue
        if '+' in l[1]:
            d = float(l[1].replace('+',''))
            if dose < d: continue
        elif '-' in l[1]:
            d = float(l[1].replace('-',''))
            if dose > d: continue
        else:
            d = float(l[1])
            if dose != d: continue
        if int(l[2]) != -999:
            low += float(l[2])
        if int(l[3]) != -999:
            high += float(l[3])
    return low, high

def readCustomRampFile(args):
    global list_CRF
    f = open(f'{args.configdir}/ramp_custom.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    list_CRF = [l.replace(' ','').split('&') for l in lines]
    for l in list_CRF:
        if len(l) < 6 or len(l) > 7:
            print(f"ERROR in file ramp_custom.txt: wrong number of arguments in entry below\n {l}")
            sys.exit()
        if l[1] != 'MOShalf' and l[1] != 'MOS2000':
            print(f"ERROR in file ramp_custom.txt: structure in entry below not recognised\n {l}")
            sys.exit()
        if l[3].lower() != 'rel' and l[3].lower() != 'abs':
            print(f"ERROR in file ramp_custom.txt: please specify rel/abs correctly in the entry below\n {l}")
            sys.exit()
        if len(l) > 6:
            if l[6].lower() != 'true' and l[6].lower() != 'false':
                print(f"ERROR in file ramp_custom.txt: please specify true/false (or nothing) for freq in the entry below\n {l}")
                sys.exit()

    return

def getCustomRamp(sample, structure, dose):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) > 6:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel', l[6].lower()=='true'
            else:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel', None
    return

def isCustomRamp(sample, structure, dose, freq):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) < 7:
                return True
            if l[6].lower() == 'true' and freq:
                return True
            if l[6].lower() == 'false' and not freq:
                return True
    return False

def getCustomRampValues(sample, structure, dose, freq):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) < 7:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
            if l[6].lower() == 'true' and freq:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
            if l[6].lower() == 'false' and not freq:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
    return None

def getGoodDirectory(sample, structure, dose, freq=False):

    if 'MOS' in structure:
        ls = good_MOS
    elif 'GCD' in structure:
        ls = good_GCD
    else:
        print(f"invalid structure: {structure}")
        return

    name = '{}_{}_{}kGy'.format(sample,structure,dose)
    if freq and 'MOS' in structure:
        name += '_1kHz'
    for l in ls:
        if not freq and 'kHz' in l:
            continue
        if name in l:
            return l
    if dose == 0:
        name = name.replace('_0kGy','_pre')
        for l in ls:
            if not freq and 'kHz' in l:
                continue
            if name in l:
                return l
    return 

def formPath(args, sample, structure, dose, freq):
    gd = getGoodDirectory(sample,structure,dose)
    if gd is None:
        return
    inpath = '{}/{}'.format(args.indir, gd)
    if 'MOS' in structure:
        meas_type = "CV"
    elif 'GCD' in structure:
        meas_type = "IV"
    else:
        print(f"invalid structure: {structure}")
        return
    inpath += '/{}_{}_needle_pad_1.txt'.format(inpath.split('/')[-1], meas_type)
    if not os.path.exists(inpath):
        print(f"ERROR! file does not exist: {inpath}")
        print("please check\n")
        return
    return inpath


def readFile(path):
    f = open(path,'r')
    l = f.read().splitlines()
    for i, line in enumerate(l):
        if '#Channel' in line: break
    l = l[i+2:]
    V = []
    IC = []
    eIC = []
    for line in l:
        line = line.split('\t')
        V.append(float(line[0]))
        IC.append(float(line[2]))
        eIC.append(float(line[3]))
    return np.array(V), np.array(IC), np.array(eIC)

def makePlot(args, sample, structure, dose, freq=False):
    path = formPath(args, sample, structure, dose, freq)
    if path is None:
        return
    V, IC, eIC = readFile(path)
    if 'GCD' in structure:
        IC = IC * -1
    V = V * -1
    return V, IC, eIC


def fitVfb(args, sample, structure, dose, Cox, freq=False):

    import numpy as np
    VIC = makePlot(args, sample, structure, dose, freq)
    if VIC is None:
        return
    V, IC, eIC = VIC
    # Remove points if needed
    n_remove = pointsToRemove(sample, structure, dose)
    if n_remove > 0:
        V = V[:-n_remove]
        IC = IC[:-n_remove]
        eIC = eIC[:-n_remove]

    # Optionally cut graph (not implemented in matplotlib, so skip for now)
    cut = True
    if sample == '3009_LR' and structure == 'MOShalf' and dose == 70 and freq == True:
        cut = False
    # TODO: implement cutGraph equivalent for matplotlib data arrays if needed

    minC = min(IC)
    maxC = max(IC)
    maxV = max(V)
    diff = maxC - minC

    # incomplete curves
    if structure == 'MOS2000' and maxC < 200:
        return

    def findX_matplotlib(yval, V, IC):
        # Find V where IC crosses yval (linear interp)
        for i in range(1, len(IC)):
            if (IC[i-1] <= yval and IC[i] >= yval) or (IC[i-1] >= yval and IC[i] <= yval):
                # linear interpolation
                x0, y0 = V[i-1], IC[i-1]
                x1, y1 = V[i], IC[i]
                if y1 == y0:
                    return x0
                return x0 + (yval-y0)*(x1-x0)/(y1-y0)
        return V[0]

    low_ramp = findX_matplotlib(minC + 0.5 * diff, V, IC)
    high_ramp = findX_matplotlib(minC + 0.9 * diff, V, IC)

    if dose == 0:
        low_ramp = findX_matplotlib(minC + 0.3 * diff, V, IC)
        high_ramp = findX_matplotlib(minC + 0.7 * diff, V, IC)

    isCustom = isCustomRamp(sample, structure, dose, freq)
    # adjusting custom ranges
    if isCustom:
        c_low, c_high, isRel = getCustomRampValues(sample, structure, dose, freq)
        if isRel:
            if c_low != -999:
                low_ramp = findX_matplotlib(minC + c_low * diff, V, IC)
            if c_high != -999:
                high_ramp = findX_matplotlib(minC + c_high * diff, V, IC)
        else:
            if c_low != -999:
                low_ramp = c_low
            if c_high != -999:
                high_ramp = c_high

    # Fit a line to the ramp region
    fit_mask = [(v >= low_ramp and v <= high_ramp) for v in V]
    fit_V = np.array([v for v, m in zip(V, fit_mask) if m])
    fit_IC = np.array([ic for ic, m in zip(IC, fit_mask) if m])
    if len(fit_V) < 2:
        return
    coeffs = np.polyfit(fit_V, fit_IC, 1)
    slope, intercept = coeffs[0], coeffs[1]

    Vfb = -1. * (intercept - Cox) / slope

    # Plotting
    nameNoExt = f"{args.outdir}/fit_{sample}_{structure}_{dose}kGy"
    plt.figure()
    plt.errorbar(V, IC, yerr=eIC, fmt='o', label='Data')
    # Set x- and y-limits with margin
    margin_x = (max(V) - min(V)) * 0.05
    margin_y = (max(IC) - min(IC)) * 0.1
    plt.xlim(min(V) - margin_x, max(V) + margin_x)
    plt.ylim(min(IC) - margin_y, max(IC) + margin_y)
    # plot fit line in the fit region
    V_fit_range = np.linspace(low_ramp, high_ramp, 100)
    ramp_fit = slope * V_fit_range + intercept
    plt.plot(V_fit_range, ramp_fit, label='Fit')
    # plot extrapolations
    V_ramp_ext = np.linspace(high_ramp, Vfb+5, 100)
    ramp_ext = slope * V_ramp_ext + intercept
    plt.plot(V_ramp_ext, ramp_ext, '--', color='blue', alpha=0.7)
    V_plat_ext = np.linspace(Vfb-10, maxV, 100)
    plat_ext = np.ones_like(V_plat_ext) * Cox
    plt.plot(V_plat_ext, plat_ext, '--', color='blue', alpha=0.7)
    # Vfb vertical line
    plt.axvline(Vfb, color='green', linestyle='--', label='Vfb')
    plt.legend(loc='best')
    plt.xlabel('-V gate [V]')
    plt.ylabel('MOS C [pF]' if 'MOS' in structure else 'diode I [nA]')
    plt.title('CV' if 'MOS' in structure else 'IV')
    if freq:
        plt.savefig(f"{nameNoExt}_1kHz.png")
    else:
        plt.savefig(f"{nameNoExt}.png")
    plt.close()
    return Vfb


def calculate_parameters(V_fb, structure, C_ox, sample):

    A = cnst.A_MOShalf
    if structure == 'MOS2000':
        A = cnst.A_MOS2000
        if C_ox > 700: 
            A = cnst.A_MOS_6inch
    elif structure == 'MOSc1' or structure == 'MOSc2':
        A = cnst.A_MOS_tracker_circle
    elif structure == 'MOS2s1' or structure == 'MOSs2' or structure == 'MOS':
        A = cnst.A_MOS_tracker_square

    if sample == '1112_LR': #cables inverted
        A = cnst.A_MOS2000
        if structure == 'MOS2000':
            A = cnst.A_MOShalf
        
    A *= 1E-6 # in m^2
    C_ox *= 1E-12 # in F

    phi_s = cnst.Chi + cnst.Eg/2. + cnst.KbTn20*ROOT.log(cnst.NA/cnst.ni)
    phi_ms = cnst.phi_m - phi_s

    N_ox = C_ox / (A*cnst.q) * (phi_ms + V_fb) # 1/m^2
    N_ox *= 1E-04 # in 1/cm^2

    t_ox = cnst.e0 * cnst.er *A / C_ox
    t_ox *= 1E09 # in nm

    return [N_ox, t_ox]


def processMOS(args, sample, structure, Cox, freq=False):

    if sample in args.MOS_exclude:
        return

    doses = []
    Vfb_values = []
    Nox_values = []
    for dose in args.doses:
        if isMOSexcluded(sample, structure, dose):
            continue
        Vfb = fitVfb(args, sample, structure, dose, Cox, freq)
        if Vfb is None:
            continue
        Nox, tox = calculate_parameters(Vfb, structure, Cox-cnst.approx_openC, sample)
        doses.append(dose)
        Vfb_values.append(Vfb)
        Nox_values.append(Nox)

    # Save to ROOT file as before (keep for compatibility)
    if len(doses) > 0:
        from array import array
        rfName = f"{args.outfiles}/dose_{sample}_{structure}"
        if freq:
            rfName += "_1kHz"
        tf = safeOpenFile(f"{rfName}.root", mode='recreate')
        gVfb = ROOT.TGraph(len(doses), array('d', doses), array('d', Vfb_values))
        gNox = ROOT.TGraph(len(doses), array('d', doses), array('d', Nox_values))
        gVfb.SetName('gVfb')
        gNox.SetName('gNox')
        gVfb.Write()
        gNox.Write()
        tf.Close()

    # Plot with matplotlib
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(doses, Vfb_values, 'o-', label='Vfb')
    plt.xlabel('dose [kGy]')
    plt.ylabel('V flat-band [-V]')
    plt.title(f'{sample} {structure}')
    plt.grid(True)
    cName = f"{args.outdir}/dose_{sample}_{structure}"
    if freq:
        plt.savefig(f"{cName}_Vfb_1kHz.png")
    else:
        plt.savefig(f"{cName}_Vfb.png")
    plt.close()

    plt.figure()
    plt.plot(doses, Nox_values, 'o-', label='Nox')
    plt.xlabel('dose [kGy]')
    plt.ylabel('oxide charge density [1/cm2]')
    plt.title(f'{sample} {structure}')
    plt.grid(True)
    if freq:
        plt.savefig(f"{cName}_Nox_1kHz.png")
    else:
        plt.savefig(f"{cName}_Nox.png")
    plt.close()

    return

def getCox(sample,structure):
    if sample in args.MOS_exclude:
        return 0
    Cox = 0
    for dose in args.doses:
        if isMOSexcluded(sample, structure, dose):
            continue
        result = makePlot(args, sample, structure, dose)
        if result is None:
            continue
        _, IC, _ = result
        maxC = max(IC)
        Cox = max(maxC, Cox)
    return Cox

def calculate_GCD_parameters(I,sample):
    I *= 1E-09 # in A
    A = cnst.A_GCD
    if sample == '3009_LR' or sample == '3010_LR':
        A = cnst.A_GCD_6inch
    elif '_SE_' in sample:
        A = cnst.A_GCD_tracker
    A *= 1E-6 # in m^2
    ni = cnst.ni * 1E+06 # in m^-3
    J = I / (cnst.q*ni*A) # in m/s
    J *= 100 # in cm/s
    return J 

def removeBadPoints(V, IC, eIC, threshold):
    import numpy as np
    # Remove points with large relative errors
    if np.max(IC) == np.min(IC):
        return V, IC, eIC
    relerr = eIC / (np.max(IC) - np.min(IC))
    mask = relerr <= threshold
    return V[mask], IC[mask], eIC[mask]

def findApproxDepletion(sample, dose, V, IC, eIC):
    import numpy as np
    baseline = np.min(IC)
    IC_base = IC - baseline
    maxY = np.max(IC_base)
    # Find first index where IC_base > 0.3 * maxY
    i = np.argmax(IC_base > maxY * 0.3)
    xL = V[i] - 10
    xH = V[i] + 10
    xL, xH = getGCDrange(sample, dose, xL, xH)
    mask = (V >= xL) & (V <= xH)
    if not np.any(mask):
        # fallback: return default values if mask is empty
        return [np.nan, np.nan, np.nan, np.nan, np.nan]
    # Find min and max IC_base in the mask region
    idx_mask = np.where(mask)[0]
    IC_base_mask = IC_base[mask]
    V_mask = V[mask]
    eIC_mask = eIC[mask]
    idx_min = np.argmin(IC_base_mask)
    idx_max = np.argmax(IC_base_mask)
    xm = V_mask[idx_min]
    ym = IC_base_mask[idx_min] + baseline
    em = eIC_mask[idx_min]
    xM = V_mask[idx_max]
    yM = IC_base_mask[idx_max] + baseline
    eM = eIC_mask[idx_max]
    return [xm, ym, xM, yM, np.sqrt(em**2 + eM**2)]

def cutGCDcurve(V, IC, eIC, cut):
    # Remove points from the end where V[-1] > cut
    while len(V) > 0 and V[-1] > cut:
        V = V[:-1]
        IC = IC[:-1]
        eIC = eIC[:-1]
    return V, IC, eIC

def cutGCDcurveLow(V, IC, eIC, cut):
    # Remove points from the start where V[0] < cut
    while len(V) > 0 and V[0] < cut:
        V = V[1:]
        IC = IC[1:]
        eIC = eIC[1:]
    return V, IC, eIC

def getGCDcurrent(args, sample, dose):
    import numpy as np
    result = makePlot(args, sample, 'GCD', dose)
    if result is None:
        return
    V, IC, eIC = result
    # Remove points if needed
    n_remove = pointsToRemove(sample, 'GCD', dose)
    if n_remove > 0:
        V = V[:-n_remove]
        IC = IC[:-n_remove]
        eIC = eIC[:-n_remove]
    low, high, gmin, gmax = getGCDcuts(sample, dose)
    # Cut low and high
    if low != -999:
        V, IC, eIC = cutGCDcurveLow(V, IC, eIC, low)
    if high != -999:
        V, IC, eIC = cutGCDcurve(V, IC, eIC, high)
    # Save original for plotting (before removeBadPoints)
    V_orig, IC_orig, eIC_orig = V[:], IC[:], eIC[:]
    # Remove bad points
    threshold = 0.05
    if sample == '3009_LR' and dose == 2:
        threshold = 0.1
    if sample == '1006_LR' and dose == 1:
        threshold = 0.1
    V, IC, eIC = removeBadPoints(V, IC, eIC, threshold)
    # Find depletion points (on cleaned data)
    xm, ym, xM, yM, e = findApproxDepletion(sample, dose, V, IC, eIC)
    # Plot using matplotlib
    import matplotlib.pyplot as plt
    plt.figure()
    # original curve (before removeBadPoints) in red
    plt.errorbar(V_orig, IC_orig, yerr=eIC_orig, fmt='o', color='red', label='Original')
    # cleaned curve in blue
    plt.errorbar(V, IC, yerr=eIC, fmt='o', color='blue', label='Cleaned')
    # depletion points as green 'o' markers (no connecting line)
    plt.errorbar([xm, xM], [ym, yM], yerr=[0,0], fmt='o', color='green', label='Depletion points')
    plt.xlabel('-V [V]')
    plt.ylabel('diode I [nA]')
    plt.title(f'{sample} GCD {dose}kGy')
    plt.legend(loc='best')
    plt.grid(True)
    plt.savefig('{}/fit_{}_GCD_{}kGy.png'.format(args.outdir, sample, dose))
    plt.close()
    return [yM - ym, e]

def processGCD(args, sample):
    print(f"Processing sample: {sample}, structure: GCD")
    if sample in args.GCD_exclude:
        return

    doses = []
    I_values = []
    J_values = []
    eI_values = []
    eJ_values = []
    for dose in args.doses:
        if dose == 0:
            continue
        if isGCDexcluded(sample, dose):
            continue
        Ie = getGCDcurrent(args, sample, dose)
        if Ie is None:
            continue
        I = Ie[0]
        e = Ie[1]
        J = calculate_GCD_parameters(I, sample)
        doses.append(dose)
        I_values.append(I)
        J_values.append(J)
        eI_values.append(e)
        eJ_values.append(e * J / I if I != 0 else 0)

    # Save to ROOT file as before (keep for compatibility)
    if len(doses) > 0:
        from array import array
        tf = safeOpenFile(f"{args.outfiles}/dose_{sample}_GCD.root", mode='recreate')
        gI = ROOT.TGraphErrors(len(doses), array('d', doses), array('d', I_values), array('d', [0]*len(doses)), array('d', eI_values))
        gJ = ROOT.TGraphErrors(len(doses), array('d', doses), array('d', J_values), array('d', [0]*len(doses)), array('d', eJ_values))
        gI.SetName('gI')
        gJ.SetName('gJ')
        gI.Write()
        gJ.Write()
        tf.Close()

    # Plot with matplotlib
    import matplotlib.pyplot as plt
    plt.figure()
    plt.errorbar(doses, I_values, yerr=eI_values, fmt='o-', label='GCD current')
    plt.xlabel('dose [kGy]')
    plt.ylabel('GCD current [nA]')
    plt.title(f'{sample} GCD')
    plt.grid(True)
    plt.savefig('{}/dose_{}_GCD_I.png'.format(args.outdir, sample))
    plt.close()

    plt.figure()
    plt.errorbar(doses, J_values, yerr=eJ_values, fmt='o-', label='surface velocity')
    plt.xlabel('dose [kGy]')
    plt.ylabel('surface velocity [cm/s]')
    plt.title(f'{sample} GCD')
    plt.grid(True)
    plt.savefig('{}/dose_{}_GCD_J.png'.format(args.outdir, sample))
    plt.close()

    return

def processSample(args, sample):
    # some customization of the structure name
    structures = ['MOShalf','MOS2000']
    if '_E_' in sample:
        structures = ['MOSc1','MOSc2','MOSs2']
    elif '_SE_' in sample:
        structures = ['MOS']

    if not args.skipMOS:
        for structure in structures:
            print(f"Processing sample: {sample}, structure: {structure}")
            if args.isAnnealing:
                processMOSannealing(args, sample, structure)
            else:
                Cox = getCox(sample,structure)
                processMOS(args, sample, structure, Cox)
                if sample == '3009_LR':
                    processMOS(args, sample, structure, Cox, freq=True)
    if not args.skipGCD and not '_E_' in sample:
        processGCD(args, sample)
    return

def initialize(args):

    readGoodList(args)
    checkGoodList(args, good_GCD)
    checkGoodList(args, good_MOS)
    readCustomRampFile(args)
    readMOSexcludeDose(args)
    readGCDexcludeDose(args)
    readRemovePoints(args)
    readGCDcuts(args)
    readGCDranges(args)

    return

def main(args):

    createPlotDirAndCopyPhp(args.outdir)
    createPlotDirAndCopyPhp(args.outfiles, copyPHP=False)

    initialize(args)

    for sample in args.samples:
        processSample(args, sample)

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",  type=str, 
                        help="Input folder")
    parser.add_argument("-o", "--outdir", type=str, 
                        help="Output folder for plots")
    parser.add_argument("-c", "--configdir", type=str, 
                        default="./configs_LabView", help="Folder with configuration files")
    parser.add_argument(      "--outfiles", type=str, 
                        help="Output folder for root files (not needed for now for annealing)")
    parser.add_argument("--doses", nargs="*", type=int,
                        help="Comma separated list of irradiation doses in kGy")
    # samples
    parser.add_argument("--samples", nargs="+", help="List of samples to be used")
    # GCD to exclude from sample list above 
    parser.add_argument("--GCD-exclude", dest="GCD_exclude", nargs="*", help="List of samples to be used")
    # MOS to exclude from sample list above
    parser.add_argument("--MOS-exclude", dest="MOS_exclude", nargs="*", help="List of samples to be used")
    parser.add_argument("--skipMOS", action="store_true", help="Skip processing MOS structures")
    parser.add_argument("--skipGCD", action="store_true", help="Skip processing GCD structures")
    parser.add_argument("-a", "--is-annealing", dest="isAnnealing", action="store_true", help="Run code on annealing data")
    parser.add_argument("--annealing-path-regexp", dest="annealingPathRegexp", type=str, default=".*", help="Use this regex to filter subfolders to be used for annealing plots")
    args = parser.parse_args()

    # Restore defaults if not given
    if args.indir is None: args.indir = indir
    if args.outdir is None: args.outdir = outdir
    if args.outfiles is None: args.outfiles = outfiles
    if args.doses is None: args.doses = doses
    if args.samples is None: args.samples = samples
    if args.GCD_exclude is None: args.GCD_exclude = GCD_exclude
    if args.MOS_exclude is None: args.MOS_exclude = MOS_exclude

    if args.isAnnealing:
        print()
        print("Warning: doing annealing, thus ignoring option --outfiles")
        print()
        args.outfiles = ""

    if args.configdir.endswith("/"):
        args.configdir = args.configdir.rstrip("/")

    main(args)



