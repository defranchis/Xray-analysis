## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *
from utility import *

from glob import glob
import os
import constants as cnst

from functionsAnnealing import *

import argparse

# now passed with options
# indir = '/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/' # now in Results_Xray/, no longer in Data/
# outdir = 'allplots'
# outfiles = 'Vfb_files'

# doses = [0,1,2,5,10,15,20,30,40,70,100,181,200,394,436,509,749,762,1030]

# samples = ['1006_LR','1008_LR','1009_LR','1010_UL','1011_LR','1003_LR','1113_LR','3009_LR',
#            '3001_UL','1112_LR','3003_UL','3103_LR','1109_LR','1105_LR','3101_LR',
#            '3010_LR','24_E_MOS','23_SE_GCD','N0538_25_LR','3007_UL','1012_UL']

# GCD_exclude = ['1008_LR','1113_LR','1105_LR','1112_LR']
# MOS_exclude = ['1012_UL']

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
    return [V, IC, eIC]

def makePlot(args, sample, structure, dose, freq=False):
    path = formPath(args, sample, structure, dose, freq)
    if path is None:
        return
    V, IC, eIC = readFile(path)
    tge = ROOT.TGraphErrors()
    for i in range(0, len(V)):
        if 'GCD' in structure:
            IC[i] *= -1.
        tge.SetPoint(i, -1*V[i], IC[i])
        tge.SetPointError(i, 0, eIC[i])
    if 'MOS' in structure:
        tge.SetName('gC')
        tge.SetTitle('CV; -V gate [V]; MOS C [pF]')
    elif 'GCD' in structure:
        tge.SetName('gI')
        tge.SetTitle('IV; -V gate [V]; diode I [nA]')
    return tge


def fitVfb(args, sample, structure, dose, Cox, freq=False):

    tge = makePlot(args, sample, structure, dose, freq)
    if tge == None: return

    for k in range(0, pointsToRemove(sample, structure, dose)):
        tge.RemovePoint(tge.GetN()-1)

    cut = True
    if sample == '3009_LR' and structure == 'MOShalf' and dose == 70 and freq == True:
        cut = False

    if cut: 
        tge = cutGraph(tge)
    
    c = ROOT.TCanvas()
    
    
    minC = min(list(tge.GetY()))
    maxC = max(list(tge.GetY()))
    maxV = max(list(tge.GetX()))
    diff = maxC-minC

    #incomplete curves
    if structure =='MOS2000' and maxC<200:
        return

    low_ramp  = findX(minC + 0.5 * diff, tge)
    high_ramp = findX(minC + 0.9 * diff, tge)

    if dose == 0:
        low_ramp  = findX(minC + 0.3 * diff, tge)
        high_ramp = findX(minC + 0.7 * diff, tge)

    isCustom = isCustomRamp(sample, structure, dose, freq)

    #adjusting custom ranges
    if isCustom:
        c_low, c_high, isRel = getCustomRampValues(sample, structure, dose, freq)
        if isRel:
            if c_low != -999:
                low_ramp = findX (minC + c_low * diff, tge)
            if c_high != -999:
                high_ramp = findX (minC + c_high * diff, tge)
        else:
            if c_low != -999:
                low_ramp = c_low
            if c_high != -999:
                high_ramp = c_high

    ramp = ROOT.TF1('ramp', 'pol1(0)', low_ramp, high_ramp)
    tge.Fit(ramp, 'q', '', low_ramp, high_ramp)

    # Vfb = (ramp.GetParameter(0)-plat.GetParameter(0))/(ramp.GetParameter(1)-plat.GetParameter(1))*(-1.)
    Vfb = -1. * (ramp.GetParameter(0)-Cox) / ramp.GetParameter(1)

    tge.Draw('ap')
    ramp.Draw('l same')

    ramp_ext = ROOT.TF1('ramp_ext', 'pol1(0)', high_ramp, Vfb+5)
    plat_ext = ROOT.TF1('plat_ext', 'pol0(0)', Vfb-10, maxV)
    ramp_ext.SetParameters(ramp.GetParameter(0), ramp.GetParameter(1))
    plat_ext.SetParameter(0, Cox)
    ramp_ext.SetLineColor(ROOT.kBlue)
    plat_ext.SetLineColor(ROOT.kBlue)

    ramp_ext.Draw('l same')
    plat_ext.Draw('l same')

    l = ROOT.TLine(Vfb, minC, Vfb, maxC)
    l.SetLineColor(ROOT.kGreen+1)
    l.Draw('l same')
    
    nameNoExt = f"{args.outdir}/fit_{sample}_{structure}_{dose}kGy"
    if freq:
        c.SaveAs(f"{nameNoExt}_1kHz.png")
    else:
        c.SaveAs(f"{nameNoExt}.png")
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

    gVfb = ROOT.TGraph()
    gNox = ROOT.TGraph()
    gVfb.SetName('gVfb')
    gNox.SetName('gNox')

    for dose in args.doses:
        if isMOSexcluded(sample, structure, dose):
            continue

        Vfb = fitVfb(args, sample, structure, dose, Cox, freq)
        if Vfb == None : continue
        Nox, tox = calculate_parameters(Vfb, structure, Cox-cnst.approx_openC, sample)
        # f.write('{} \t {} \t {} \n'.format(dose,Vfb,Nox))
        gVfb.SetPoint(gVfb.GetN(), dose, Vfb)
        gNox.SetPoint(gNox.GetN(), dose, Nox)
    # f.close()

    rfName = f"{args.outfiles}/dose_{sample}_{structure}"
    if freq:
        rfName += "_1kHz"
    tf = safeOpenFile(f"{rfName}.root", mode='recreate')
    gVfb.Write()
    gNox.Write()
    tf.Close()

    c = ROOT.TCanvas()
    gVfb.SetTitle('{} {}; dose [kGy]; V flat-band [-V]'.format(sample,structure))
    gVfb.SetMarkerStyle(7)
    gVfb.Draw('apl')
    cName = f"{args.outdir}/dose_{sample}_{structure}"
    if freq:
        c.SaveAs(f"{cName}_Vfb_1kHz.png")
    else:
        c.SaveAs(f"{cName}_Vfb.png")
    c.Clear()
    gNox.SetTitle('{} {}; dose [kGy]; oxide charge density [1/cm2]'.format(sample,structure))
    gNox.SetMarkerStyle(7)
    gNox.Draw('apl')
    if freq:
        c.SaveAs(f"{cName}_Nox_1kHz.png")
    else:
        c.SaveAs(f"{cName}_Nox.png")

    return

def getCox(sample,structure):
    if sample in args.MOS_exclude:
        return 0
    Cox = 0
    for dose in args.doses:
        if isMOSexcluded(sample, structure, dose):
            continue
        tge = makePlot(args, sample, structure, dose)
        if tge is None:
            continue
        maxC = max(list(tge.GetY()))
        Cox = max(maxC,Cox)
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

def removeBadPoints(tge,threshold):
    eY = list(tge.GetEY())
    Y = list(tge.GetY())
    # relerr = [i/j for  i,j in zip(eY,Y)]
    if max(Y)==min(Y):
        return tge
    relerr = [i/(max(Y)-min(Y)) for  i in eY]

    for i, re in enumerate(relerr):
        if re > threshold:
            tge.RemovePoint(i)
            break
    if i < len(relerr)-1:
        tge = removeBadPoints(tge,threshold)
    return tge

def findApproxDepletion(sample,dose,tge):
    X = list(tge.GetX())
    Y = list(tge.GetY())
    eY = list(tge.GetEY())
    baseline = min(Y)
    Y = [y-baseline for y in Y]
    for i,y in enumerate(Y):
        if y > max(Y)*0.3: 
            break
    
    xL = X[i]-10
    xH = X[i]+10
    xL, xH = getGCDrange(sample,dose,xL,xH)

    ym = +999
    yM = -999
    xm = -999
    xM = +999
    eM = 0
    em = 0

    for i,x in enumerate(X):
        if x < xL or x > xH:
            continue
        if Y[i] > yM: 
            yM = Y[i]
            eM = eY[i]
            xM = x
        if Y[i] < ym: 
            ym = Y[i]
            em = eY[i]
            xm = x
    return [xm,ym+baseline,xM,yM+baseline,(em**2+eM**2)**0.5]

def cutGCDcurve(tge,cut):
    if list(tge.GetX())[-1] > cut:
        tge.RemovePoint(tge.GetN()-1)
    if list(tge.GetX())[-1] > cut:
        tge = cutGCDcurve(tge,cut)

    tge.SetMaximum(max(list(tge.GetY()))*1.1)
    return tge

def cutGCDcurveLow(tge,cut):
    if list(tge.GetX())[0] < cut:
        tge.RemovePoint(0)
    if list(tge.GetX())[0] < cut:
        tge = cutGCDcurveLow(tge,cut)

    tge.SetMaximum(max(list(tge.GetY()))*1.1)
    return tge

def getGCDcurrent(args, sample, dose):
    tge = makePlot(args, sample, 'GCD', dose)
    if tge == None: return

    for k in range(0,pointsToRemove(sample,'GCD',dose)):
        tge.RemovePoint(tge.GetN()-1)
        
    low, high, gmin, gmax = getGCDcuts(sample,dose)

    if low != -999:
        tge = cutGCDcurveLow(tge,low)
    if high != -999:
        tge = cutGCDcurve(tge,high)
    if gmin != -999:
        tge.SetMinimum(gmin)
    if gmax != -999:
        tge.SetMaximum(gmax)

    tge_orig = tge.Clone()    
    threshold = 0.05
    if sample == '3009_LR' and dose == 2:
        threshold = 0.1
    if sample == '1006_LR' and dose == 1:
        threshold = 0.1
    tge = removeBadPoints(tge,threshold)

    xm, ym, xM, yM, e = findApproxDepletion(sample,dose,tge)

    at = ROOT.TGraph()
    at.SetPoint(0,xm,ym)
    at.SetPoint(1,xM,yM)
    at.SetMarkerStyle(8)
    at.SetMarkerColor(ROOT.kGreen)
        
    c = ROOT.TCanvas()
    tge_orig.SetLineColor(ROOT.kRed)
    tge_orig.SetMarkerColor(ROOT.kRed)
    tge_orig.Draw('apl')
    tge.Draw('pl same')
    at.Draw('p same')
    c.SaveAs('{}/fit_{}_GCD_{}kGy.png'.format(args.outdir,sample,dose))

    return [yM-ym,e]

def processGCD(args, sample):
    
    if sample in args.GCD_exclude:
        return

    gI = ROOT.TGraphErrors()
    gJ = ROOT.TGraphErrors()
    gI.SetName('gI')
    gJ.SetName('gJ')

    for dose in args.doses:
        if dose==0: continue
        if isGCDexcluded(sample,dose):
            continue
        
        Ie = getGCDcurrent(args, sample, dose)
        if Ie == None : continue
        I = Ie[0]
        e = Ie[1]
        J = calculate_GCD_parameters(I,sample)
        gI.SetPoint(gI.GetN(),dose,I)
        gJ.SetPoint(gJ.GetN(),dose,J)
        gI.SetPointError(gI.GetN()-1, 0, e)
        gJ.SetPointError(gJ.GetN()-1, 0, e*J/I)


    tf = safeOpenFile(f"{args.outfiles}/dose_{sample}_GCD.root", mode='recreate')
    gI.Write()
    gJ.Write()
    tf.Close()

    c = ROOT.TCanvas()
    gI.SetTitle('{} GCD; dose [kGy]; GCD current [nA]'.format(sample))
    gI.SetMarkerStyle(7)
    gI.Draw('apl')
    c.SaveAs('{}/dose_{}_GCD_I.png'.format(args.outdir, sample))
    c.Clear()
    gJ.SetTitle('{} GCD; dose [kGy]; surface velocity [cm/s]'.format(sample))
    gJ.SetMarkerStyle(7)
    gJ.Draw('apl')
    c.SaveAs('{}/dose_{}_GCD_J.png'.format(args.outdir,sample))

    return

def processSample(args, sample):
    # some customization of the structure name
    structures = ['MOShalf','MOS2000']
    if '_E_' in sample:
        structures = ['MOSc1','MOSc2','MOSs2']
    elif '_SE_' in sample:
        structures = ['MOS']

    if not args.skipStructure == "MOS":
        for structure in structures:
            if args.isAnnealing:
                processMOSannealing(args, sample, structure)
            else:
                Cox = getCox(sample,structure)
                processMOS(args, sample, structure, Cox)
                if sample == '3009_LR':
                    processMOS(args, sample, structure, Cox, freq=True)
    if not args.skipStructure == "GCD" and not '_E_' in sample:
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
                        default="/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/", help="Input folder")
    parser.add_argument("-o", "--outdir", type=str, 
                        default="./allplots", help="Output folder for plots")
    parser.add_argument("-c", "--configdir", type=str, 
                        default="./configs", help="Folder with configuration files")
    parser.add_argument(      "--outfiles", type=str, 
                        default="./Vfb_files", help="Output folder for root files (not needed for now for annealing)")
    parser.add_argument("--doses", nargs="*", type=int, default=[0,1,2,5,10,15,20,30,40,70,100,181,200,394,436,509,749,762,1030],
                        help="Comma separated list of irradiation doses in kGy")
    # samples
    parser.add_argument("--samples", nargs="+", default=['1006_LR','1008_LR','1009_LR','1010_UL','1011_LR','1003_LR','1113_LR','3009_LR','3001_UL','1112_LR','3003_UL','3103_LR','1109_LR','1105_LR','3101_LR','3010_LR','24_E_MOS','23_SE_GCD','N0538_25_LR','3007_UL','1012_UL'], help="List of samples to be used")
    # GCD to exclude from sample list above 
    # FIXME: nargs="+" will request at least one argument, so the default should be empty array if one does not to use the options, will change it once we have these input options stored somewhere (e.g. passed from external configuration files)
    parser.add_argument("--GCD-exclude", dest="GCD_exclude", nargs="*", default=['1008_LR','1113_LR','1105_LR','1112_LR'], help="List of samples to be used")
    # MOS to exclude from sample list above
    # FIXME: see previous comment
    parser.add_argument("--MOS-exclude", dest="MOS_exclude", nargs="*", default=['1012_UL'], help="List of samples to be used")
    parser.add_argument("--skip", dest="skipStructure", choices=["", "MOS", "GCD"],
                        default="", help="To exclude completely one kind of structure and be faster")
    parser.add_argument("-a", "--is-annealing", dest="isAnnealing", action="store_true", help="Run code on annealing data")
    parser.add_argument("--annealing-path-regexp", dest="annealingPathRegexp", type=str, default=".*", help="Use this regex to filter subfolders to be used for annealing plots")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    # print(args.doses)
    # for x in args.doses:
    #     if not isinstance(x, int):
    #         print("%s is not integer!" % str(x))
    # print(args.samples)
    # quit()

    # print(args.GCD_exclude)
    # print(f"len = {len(args.GCD_exclude)}")
    # quit()

    if args.isAnnealing:
        print()
        print("Warning: doing annealing, thus ignoring option --outfiles")
        print()
        args.outfiles = ""

    if args.configdir.endswith("/"):
        args.configdir = args.configdir.rstrip("/")

    main(args)



