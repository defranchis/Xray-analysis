import ROOT
from ROOT import *
from fitVfb import cutGraph, removeBadPoints, findApproxDepletion, findX, cutGCDcurveLow
import os
from glob import glob

gROOT.SetBatch(True)

outdir = 'allplots'
if not os.path.exists(outdir):                                            
    os.makedirs(outdir) 

outfiles = 'Vfb_files'
if not os.path.exists(outfiles):                                            
    os.makedirs(outfiles) 

doses = [0,1,2,5,10,20]

def getPath(sample,structure,dose):
    path = '/eos/user/h/hgsensor/HGCAL_test_results/Data/logs_obelix_setup/{}_m20C_{}kGy'.format(sample,dose)
    path = glob('{}/*/'.format(path))[-1]
    if 'MOS' in structure:
        path += 'cv_{}_m20C_{}kGy_{}.dat'.format(sample,dose,structure)
    elif 'GCD' in structure:
        path += 'iv_{}_m20C_{}kGy_{}.dat'.format(sample,dose,structure)
    if not os.path.exists(path):
        print path
    return path


def getPlot(sample,structure,dose):
    path = getPath(sample,structure,dose)
    f = open(path,'r')
    lines = f.read().splitlines()
    V = []
    meas = []
    err = []
    for line in lines:
        if '#' in line: continue
        data = line.split()
        V.append(-1.*float(data[0]))
        if 'MOS' in structure:
            meas.append(float(data[-3]))
            err.append(0)

        elif 'GCD' in structure:
            meas.append(-1.e09*float(data[2]))
            err.append(1e09*float(data[3]))

    tge = TGraphErrors()
    for i, v in enumerate(V):
        tge.SetPoint(i,v,meas[i])
        tge.SetPointError(i,0,err[i])
    if tge.GetN()<2: return
    return tge

def fitVfb(sample,structure,dose):
    tge = getPlot(sample,structure,dose)
    if tge == None: return

    tge = cutGraph(tge)
    
    c = TCanvas()

    minC = min(list(tge.GetY()))
    maxC = max(list(tge.GetY()))
    diff = maxC-minC

    if structure =='MOS2000' and maxC*1e12<200:
        return

    low_ramp = findX (minC + .5 * diff, tge)
    high_ramp = findX (minC + .9 * diff, tge)

    if dose == 0:
        low_ramp = findX (minC + .3 * diff, tge)
        high_ramp = findX (minC + .7 * diff, tge)


    ramp = TF1('ramp','pol1(0)',low_ramp,high_ramp)
    tge.Fit(ramp,'q','',low_ramp,high_ramp)

    low_plat = high_ramp*1.1
    high_plat = high_ramp*1.2

    if dose == 0:
        low_plat = high_ramp*3
        high_plat = high_ramp*4
        maxV = max(list(tge.GetX()))
        if low_plat > maxV:
            high_plat = maxV
            low_plat = .9*maxV

    
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
    
    c.SaveAs('{}/fit_{}_{}_{}kGy.png'.format(outdir,sample,structure,dose))
    return Vfb



def processMOS(sample,structure):
    gVfb = TGraph()
    gVfb.SetName('gVfb')

    for dose in doses:
        Vfb = fitVfb(sample,structure,dose)
        if Vfb == None : continue
        gVfb.SetPoint(gVfb.GetN(),dose,Vfb)
    
    tf = TFile.Open('{}/dose_{}_{}.root'.format(outfiles,sample,structure),'recreate')
    gVfb.Write()

    c = TCanvas()
    gVfb.SetTitle('{} {}; dose [kGy]; V flat-band [-V]'.format(sample,structure))
    gVfb.SetMarkerStyle(7)
    gVfb.Draw('apl')
    c.SaveAs('{}/dose_{}_{}_Vfb.png'.format(outdir,sample,structure))

    return

def getGCDcurrent(sample,dose):
    tge = getPlot(sample,'GCD',dose)
    if tge == None: return

    tge_orig = tge.Clone()    
    threshold = 0.15
    tge = removeBadPoints(tge,threshold)
    if sample == '3008_UL' and dose == 1:
        tge = cutGCDcurveLow(tge,-9)

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
    gI = TGraphErrors()
    gI.SetName('gI')

    for dose in doses:
        if sample == '3008_UL':
            if dose == 0 or dose >=5:
                continue
        Ie = getGCDcurrent(sample,dose)
        if Ie == None : continue
        I = Ie[0]
        e = Ie[1]
        gI.SetPoint(gI.GetN(),dose,I)
        gI.SetPointError(gI.GetN()-1,0,e)

    tf = TFile.Open('{}/dose_{}_GCD.root'.format(outfiles,sample),'recreate')
    gI.Write()

    c = TCanvas()
    gI.SetTitle('{} GCD; dose [kGy]; GCD current [nA]'.format(sample))
    gI.SetMarkerStyle(7)
    gI.Draw('apl')
    c.SaveAs('{}/dose_{}_GCD_I.png'.format(outdir,sample))

    return

def processSample(sample):
    structures = ['MOShalf','MOS2000']
    for structure in structures:
        processMOS(sample,structure)
    processGCD(sample)
    return

processSample('3008_UL')
