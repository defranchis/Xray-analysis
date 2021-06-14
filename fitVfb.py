import ROOT
from ROOT import *
from glob import glob
import os
import constants as cnst
import sys

gROOT.SetBatch(True)

indir = '/eos/user/h/hgsensor/HGCAL_test_results/Data/'
outdir = 'allplots'
outfiles = 'Vfb_files'

doses = [0,1,2,5,10,15,20,30,40,70,100,181,200,394,436,509,749,762,1030]

samples = ['1006_LR','1008_LR','1009_LR','1010_UL','1011_LR','1003_LR','1113_LR','3009_LR',
           '3001_UL','1112_LR','3003_UL','3103_LR','1109_LR','1105_LR','3101_LR',
           '3010_LR','24_E_MOS','23_SE_GCD','N0538_25_LR','3007_UL','1012_UL']

GCD_exclude = ['1008_LR','1113_LR','1105_LR']
MOS_exclude = ['1012_UL']

def checkGoodList(ls):
    for i,l1 in enumerate(ls):
        for j,l2 in enumerate(ls):
            if i==j: continue
            if l1 in l2: 
                if l1+'_1kHz' in l2: continue
                if 'GCD' in l1 and l1.split('_')[1]+'_'+l1.split('_')[2] in GCD_exclude:
                    continue
                print '\nWARNING: potential duplicate:', l1, l2
                print 'please check\n'
    return

def readGoodList():
    global good_GCD
    global good_MOS
    f_GCD = open('good_GCD.txt','r') 
    f_MOS = open('good_MOS.txt','r') 
    good_GCD = f_GCD.read().splitlines()
    good_MOS = f_MOS.read().splitlines()
    return

def readMOSexcludeDose():
    global MOS_exclude_dose
    f = open('MOS_exclude_dose.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    MOS_exclude_dose = [l.replace(' ','').split('&') for l in lines]
    for l in MOS_exclude_dose:
        if len(l) != 3:
            print 'ERROR in file MOS_exclude_dose.txt: wrong number of arguments in entry below\n {}'.format(l)
            sys.exit()
        if l[1] != 'MOShalf' and l[1] != 'MOS2000':
            print 'ERROR in file MOS_exclude_dose.txt: structure in entry below not recognised\n {}'.format(l)
            sys.exit()
    return

def isMOSexcluded(sample,structure,dose):
    for l in MOS_exclude_dose:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            return True
    return False

def readCustomRampFile():
    global list_CRF
    f = open('ramp_custom.txt','r')
    lines = f.read().splitlines()
    f.close()
    lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
    list_CRF = [l.replace(' ','').split('&') for l in lines]
    for l in list_CRF:
        if len(l) < 6 or len(l) > 7:
            print 'ERROR in file ramp_custom.txt: wrong number of arguments in entry below\n {}'.format(l)
            sys.exit()
        if l[1] != 'MOShalf' and l[1] != 'MOS2000':
            print 'ERROR in file ramp_custom.txt: structure in entry below not recognised\n {}'.format(l)
            sys.exit()
        if l[3].lower() != 'rel' and l[3].lower() != 'abs':
            print 'ERROR in file ramp_custom.txt: please specify rel/abs correctly in the entry below\n {}'.format(l)
            sys.exit()
        if len(l) > 6:
            if l[6].lower() != 'true' and l[6].lower() != 'false':
                print 'ERROR in file ramp_custom.txt: please specify true/false (or nothing) for freq in the entry below\n {}'.format(l)
                sys.exit()

    return

def getCustomRamp(sample,structure,dose):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) > 6:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel', l[6].lower()=='true'
            else:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel', None
    return

def isCustomRamp(sample,structure,dose,freq):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) < 7:
                return True
            if l[6].lower() == 'true' and freq:
                return True
            if l[6].lower() == 'false' and not freq:
                return True
    return False

def getCustomRampValues(sample,structure,dose,freq):
    for l in list_CRF:
        if l[0] == sample and l[1] == structure and float(l[2])==dose:
            if len(l) < 7:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
            if l[6].lower() == 'true' and freq:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
            if l[6].lower() == 'false' and not freq:
                return float(l[4]), float(l[5]), l[3].lower() == 'rel'
    return None

def getGoodDirectory(sample, structure, dose,freq=False):

    if 'MOS' in structure:
        ls = good_MOS
    elif 'GCD' in structure:
        ls = good_GCD
    else:
        print 'invalid structure: ', structure
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

def formPath(sample, structure, dose,freq):
    gd = getGoodDirectory(sample,structure,dose)
    if gd is None:
        return
    inpath = '{}/{}'.format(indir,gd)
    if 'MOS' in structure:
        inpath += '/{}_CV_needle_pad_1.txt'.format(inpath.split('/')[-1])
    elif 'GCD' in structure:
        inpath += '/{}_IV_needle_pad_1.txt'.format(inpath.split('/')[-1])
    else:
        print 'invalid structure: ', structure
        return
    if not os.path.exists(inpath):
        print 'ERROR! file does not exist:', inpath
        print 'please check\n'
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

def makePlot(sample,structure,dose,freq=False):
    path = formPath(sample, structure, dose,freq)
    if path is None:
        return
    V, IC, eIC = readFile(path)
    tge = TGraphErrors()
    for i in range(0,len(V)):
        if 'GCD' in structure:
            IC[i] *= -1.
        tge.SetPoint(i,-1*V[i],IC[i])
        tge.SetPointError(i,0,eIC[i])
    if 'MOS' in structure:
        tge.SetName('gC')
        tge.SetTitle('CV; -V gate [V]; MOS C [pF]')
    elif 'GCD' in structure:
        tge.SetName('gI')
        tge.SetTitle('IV; -V gate [V]; diode I [nA]')

    return tge

def graph_derivative(g):
    der = TGraph()
    for i in range(0,g.GetN()-1):
        der.SetPoint(i,(g.GetPointX(i+1)+g.GetPointX(i))/2.,(g.GetPointY(i+1)-g.GetPointY(i))/(g.GetPointX(i+1)-g.GetPointX(i)))
    return der

def cutGraph(old):
    new = TGraphErrors()
    X = list(old.GetX())
    Y = list(old.GetY())
    EX = list(old.GetEX())
    EY = list(old.GetEY())
    i = Y.index(min(Y))
    X = X[i:]
    Y = Y[i:]
    EX = EX[i:]
    EY = EY[i:]

    for i,x in enumerate(X):
        new.SetPoint(i,x,Y[i])
        new.SetPointError(i,EX[i],EY[i])

    return new


def findX(yy,tge):
    X = list(tge.GetX())
    Y = list(tge.GetY())
    for i,y in enumerate(Y):
        if y>yy: break
    return (X[i]+X[i-1])*.5


def fitVfb(sample,structure,dose,freq=False):

    tge = makePlot(sample,structure,dose,freq)
    if tge == None: return

    if sample == '1112_LR' and structure == 'MOShalf':
        tge.RemovePoint(tge.GetN()-1)
        tge.RemovePoint(tge.GetN()-1)
        tge.RemovePoint(tge.GetN()-1)

    cut = True
    if sample == '3009_LR' and structure == 'MOShalf' and dose == 70 and freq == True:
        cut = False

    if cut: 
        tge = cutGraph(tge)
    
    c = TCanvas()

    minC = min(list(tge.GetY()))
    maxC = max(list(tge.GetY()))
    diff = maxC-minC

    #incomplete curves
    if structure =='MOS2000' and maxC<200:
        return

    low_ramp = findX (minC + .5 * diff, tge)
    high_ramp = findX (minC + .9 * diff, tge)

    if dose == 0:
        low_ramp = findX (minC + .3 * diff, tge)
        high_ramp = findX (minC + .7 * diff, tge)

    isCustom = isCustomRamp(sample,structure,dose,freq)

    #adjusting custom ranges
    if isCustom:
        c_low, c_high, isRel = getCustomRampValues(sample,structure,dose,freq)
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
            

    ramp = TF1('ramp','pol1(0)',low_ramp,high_ramp)
    tge.Fit(ramp,'q','',low_ramp,high_ramp)

    low_plat = high_ramp*1.1
    high_plat = high_ramp*1.2

    if dose == 0:
        low_plat = high_ramp*1.5
        high_plat = high_ramp*1.7
        maxV = max(list(tge.GetX()))
        if low_plat > maxV:
            high_plat = maxV
            low_plat = .9*maxV
        if sample == '1008_LR' or sample == '1010_UL' or sample == '3001_UL':
            if structure == 'MOShalf':
                low_plat = .8 * maxV
        if sample == '23_SE_GCD' or sample == 'N0538_25_LR' or sample == '3007_UL':
            low_plat = 5
            high_plat = 7
    if sample == '1006_LR' and structure == 'MOS2000' and dose >=10:
        high_plat = max(list(tge.GetX()))
        low_plat = .95 * high_plat
    if sample == 'N0538_25_LR' and structure == 'MOS2000' and dose ==5:
        high_plat = max(list(tge.GetX()))
        low_plat = .95 * high_plat

    if sample == '3010_LR' and structure == 'MOS2000' and dose == 2:
        low_plat = 145
        high_ramp = 150

    if sample == '24_E_MOS' and structure == 'MOSc2' and dose == 2:
        low_plat = 145
        high_ramp = 150


    if '1009' in sample and dose == 40 and structure == 'MOShalf':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '1011_LR' in sample and dose == 20 and structure == 'MOShalf':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '3001' in sample and dose == 70 and structure == 'MOShalf':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '1112_LR' in sample and dose == 100 and structure == 'MOShalf':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '1011_LR' in sample and dose == 10 and structure == 'MOS2000':
        low_plat = high_ramp*1.03
        high_plat = high_ramp*1.1
    if '1011_LR' in sample and dose == 2 and structure == 'MOS2000':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '1008_LR' in sample and dose == 2 and structure == 'MOS2000':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '1113_LR' in sample and dose == 5 and structure == 'MOS2000':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if '3001_UL' in sample and dose == 10 and structure == 'MOS2000':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
    if sample == '1010_UL' and (dose == 5 or dose == 10) and structure == 'MOS2000':
        low_plat = high_ramp*1.05
        high_plat = high_ramp*1.1
 
    if sample == '23_SE_GCD' and dose == 70:
        high_plat = max(list(tge.GetX()))
        low_plat = high_plat *.95
    if sample == 'N0538_25_LR' and dose == 20:
        high_plat = max(list(tge.GetX()))
        low_plat = high_plat *.95
    


    plat = TF1('plat','pol1(0)',low_plat,high_plat)
    tge.Fit(plat,'q','',low_plat,high_plat)

    Vfb = (ramp.GetParameter(0)-plat.GetParameter(0))/(ramp.GetParameter(1)-plat.GetParameter(1))*(-1.)

    tge.Draw('ap')
    plat.Draw('l same')
    ramp.Draw('l same')

    ramp_ext = TF1('ramp_ext','pol1(0)',high_ramp,Vfb+5)
    plat_ext = TF1('plat_ext','pol1(0)',Vfb-10,low_plat)
    ramp_ext.SetParameters(ramp.GetParameter(0),ramp.GetParameter(1))
    plat_ext.SetParameters(plat.GetParameter(0),plat.GetParameter(1))
    ramp_ext.SetLineColor(kBlue)
    plat_ext.SetLineColor(kBlue)

    ramp_ext.Draw('l same')
    plat_ext.Draw('l same')

    l = TLine(Vfb,minC,Vfb,maxC)
    l.SetLineColor(kGreen+1)
    l.Draw('l same')
    
    if freq:
        c.SaveAs('{}/fit_{}_{}_{}kGy_1kHz.png'.format(outdir,sample,structure,dose))
    else:
        c.SaveAs('{}/fit_{}_{}_{}kGy.png'.format(outdir,sample,structure,dose))
    return Vfb


def calculate_parameters(V_fb, structure, C_ox, sample):

    A = cnst.A_MOShalf
    if structure == 'MOS2000':
        A = cnst.A_MOS2000
        if C_ox > 700: A = cnst.A_MOS_6inch
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

    phi_s = cnst.Chi + cnst.Eg/2. + cnst.KbTn20*log(cnst.NA/cnst.ni)
    phi_ms = cnst.phi_m - phi_s

    N_ox = C_ox / (A*cnst.q) * (phi_ms + V_fb) # 1/m^2
    N_ox *= 1E-04 # in 1/cm^2

    t_ox = cnst.e0 * cnst.er *A / C_ox
    t_ox *= 1E09 # in nm

    return [N_ox, t_ox]


def processMOS(sample,structure,Cox,freq=False):

    if sample in MOS_exclude:
        return

    gVfb = TGraph()
    gNox = TGraph()
    gVfb.SetName('gVfb')
    gNox.SetName('gNox')

    for dose in doses:
        if isMOSexcluded(sample,structure,dose):
            continue

        Vfb = fitVfb(sample,structure,dose,freq)
        if Vfb == None : continue
        Nox, tox = calculate_parameters(Vfb, structure, Cox, sample)
        # f.write('{} \t {} \t {} \n'.format(dose,Vfb,Nox))
        gVfb.SetPoint(gVfb.GetN(),dose,Vfb)
        gNox.SetPoint(gNox.GetN(),dose,Nox)
    # f.close()

    if freq:
        tf = TFile.Open('{}/dose_{}_{}_1kHz.root'.format(outfiles,sample,structure),'recreate')
    else:
        tf = TFile.Open('{}/dose_{}_{}.root'.format(outfiles,sample,structure),'recreate')
    gVfb.Write()
    gNox.Write()

    c = TCanvas()
    gVfb.SetTitle('{} {}; dose [kGy]; V flat-band [-V]'.format(sample,structure))
    gVfb.SetMarkerStyle(7)
    gVfb.Draw('apl')
    if freq:
        c.SaveAs('{}/dose_{}_{}_Vfb_1kHz.png'.format(outdir,sample,structure))
    else:
        c.SaveAs('{}/dose_{}_{}_Vfb.png'.format(outdir,sample,structure))
    c.Clear()
    gNox.SetTitle('{} {}; dose [kGy]; oxide charge density [1/cm2]'.format(sample,structure))
    gNox.SetMarkerStyle(7)
    gNox.Draw('apl')
    if freq:
        c.SaveAs('{}/dose_{}_{}_Nox_1kHz.png'.format(outdir,sample,structure))
    else:
        c.SaveAs('{}/dose_{}_{}_Nox.png'.format(outdir,sample,structure))

    return

def getCox(sample,structure):
    if sample in MOS_exclude:
        return 0
    tge = makePlot(sample,structure,0)
    Cox = max(list(tge.GetY()))
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
        if y > max(Y)*.3: break
    
    xL = X[i]-10
    xH = X[i]+10

    if sample == '3103_LR' and dose == 40:
        xH += 10
    if sample == '3101_LR':
        xL -= 20
        if dose > 1000:
            xL -= 10
    if sample == '1010_UL' and dose == 1:
        xL += 5
    if sample == '3009_LR' and dose >= 20:
        xL -= 15
        xH += 20
    if sample == '1006_LR' and dose >= 70:
        xL -= 20

    if sample == '3010_LR' and dose >= 40:
        xL -= 20
        xH += 20

    if sample == '3010_LR' and dose == 20:
        xH +=10
        
    if sample == '23_SE_GCD' and dose >=40:
        xH +=30
        xL -=20
    if sample == 'N0538_25_LR' and dose == 1:
        xH -= 7
        xL -= 10
    if sample == 'N0538_25_LR' and dose >= 40:
        xL -= 20
        xH += 20
    if sample == '3007_UL' and dose >= 20:
        xL -= 20

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
    return [xm,ym+baseline,xM,yM+baseline,(em**2+eM**2)**.5]

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

def getGCDcurrent(sample,dose):
    tge = makePlot(sample,'GCD',dose)
    if tge == None: return

    if sample == '1112_LR':
        tge = cutGCDcurve(tge,50)
    if sample == '3103_LR' and dose == 20:
        tge = cutGCDcurve(tge,40)
    if sample == '3103_LR' and dose == 40:
        tge = cutGCDcurve(tge,47)
    if sample == '3003_UL' and dose == 5:
        tge.RemovePoint(tge.GetN()-1)
        tge = cutGCDcurve(tge,40)
    if sample == '1012_UL' and dose == 1:
        tge.RemovePoint(tge.GetN()-1)
        tge.SetMinimum(0)
        tge.SetMaximum(.4)
    if sample == '1011_LR' and dose <= 5:
        tge = cutGCDcurve(tge,40)
    if sample == '1011_LR' and dose >= 10:
        tge = cutGCDcurve(tge,55)
    if sample == '1109_LR' and dose == 40:
        tge = cutGCDcurveLow(tge,50)
    if sample == '1109_LR' and dose == 1:
        tge = cutGCDcurveLow(tge,-2)
    if sample == '1109_LR' and dose >= 40:
        tge = cutGCDcurve(tge,80)
    if sample == '3007_UL' and dose == 20:
        tge = cutGCDcurve(tge,45)
    if sample == '3007_UL' and dose == 40:
        tge = cutGCDcurve(tge,50)



    tge_orig = tge.Clone()    
    threshold = 0.05
    if sample == '3009_LR' and dose == 2:
        threshold = 0.1
    if sample == '1006_LR' and dose == 1:
        threshold = 0.1
    tge = removeBadPoints(tge,threshold)

    xm, ym, xM, yM, e = findApproxDepletion(sample,dose,tge)

    at = TGraph()
    at.SetPoint(0,xm,ym)
    at.SetPoint(1,xM,yM)
    at.SetMarkerStyle(8)
    at.SetMarkerColor(kGreen)
        
    c = TCanvas()
    tge_orig.SetLineColor(kRed)
    tge_orig.SetMarkerColor(kRed)
    tge_orig.Draw('apl')
    tge.Draw('pl same')
    at.Draw('p same')
    c.SaveAs('{}/fit_{}_GCD_{}kGy.png'.format(outdir,sample,dose))

    return [yM-ym,e]

def processGCD(sample):
    
    if sample in GCD_exclude:
        return

    gI = TGraphErrors()
    gJ = TGraphErrors()
    gI.SetName('gI')
    gJ.SetName('gJ')

    for dose in doses:
        if dose==0: continue
        if sample == '3103_LR' and dose >=70: continue
        if sample == '3101_LR' and dose ==100: continue
        if sample == '1112_LR' and dose !=1 : continue
        if sample == '23_SE_GCD' and dose <=1: continue
        if sample == '3007_UL' and dose >=70: continue

        Ie = getGCDcurrent(sample,dose)
        if Ie == None : continue
        I = Ie[0]
        e = Ie[1]
        J = calculate_GCD_parameters(I,sample)
        gI.SetPoint(gI.GetN(),dose,I)
        gJ.SetPoint(gJ.GetN(),dose,J)
        gI.SetPointError(gI.GetN()-1,0,e)
        gJ.SetPointError(gJ.GetN()-1,0,e*J/I)


    tf = TFile.Open('{}/dose_{}_GCD.root'.format(outfiles,sample),'recreate')
    gI.Write()
    gJ.Write()

    c = TCanvas()
    gI.SetTitle('{} GCD; dose [kGy]; GCD current [nA]'.format(sample))
    gI.SetMarkerStyle(7)
    gI.Draw('apl')
    c.SaveAs('{}/dose_{}_GCD_I.png'.format(outdir,sample))
    c.Clear()
    gJ.SetTitle('{} GCD; dose [kGy]; surface velocity [cm/s]'.format(sample))
    gJ.SetMarkerStyle(7)
    gJ.Draw('apl')
    c.SaveAs('{}/dose_{}_GCD_J.png'.format(outdir,sample))

    return

def processSample(sample):
    structures = ['MOShalf','MOS2000']
    if '_E_' in sample:
        structures = ['MOSc1','MOSc2','MOSs2']
    elif '_SE_' in sample:
        structures = ['MOS']

    for structure in structures:
        Cox = getCox(sample,structure)
        processMOS(sample,structure,Cox-cnst.approx_openC)
        if sample == '3009_LR':
            processMOS(sample,structure,Cox-cnst.approx_openC,freq=True)
    if not '_E_' in sample:
            processGCD(sample)
    return


def main():

    if not os.path.exists(outdir):                                            
        os.makedirs(outdir) 
    if not os.path.exists(outfiles):                                            
        os.makedirs(outfiles) 

    readGoodList()
    checkGoodList(good_GCD)
    checkGoodList(good_MOS)
    readCustomRampFile()
    readMOSexcludeDose()

    for sample in samples:
        processSample(sample)

    return

if __name__ == "__main__":
    main()



