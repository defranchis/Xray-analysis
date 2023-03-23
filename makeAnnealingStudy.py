# example
# python3 makeAnnealingStudy.py -o plots/annealing/ -c "configs_Nov2021/" -t "temperatureLogs/temperatureLog-2022-02-25-17-25-55.dat" --max-dose 0 -s "N4789-12_UL" --last-step 52
# python3 makeAnnealingStudy.py -o plots/annealing/ -c "configs_Nov2021/" -t "temperatureLogs/temperatureLog-2022-03-04-18-39-10.dat" --max-dose 0 -s "N4789-12_LR"

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
import math
import constants as cnst
import sampleStyle as sast
import datetime

from functionsAnnealing import *

from makeXrayStudy import siliconSensorSample
from makeXrayStudy import getOxideChargeDensity

import argparse

tgeTemperature = ROOT.TGraphErrors()
tgeTemperature.SetName(f"temperature")
def temperatureGraphForFit(x, p):
    #return gr.Eval(x[0], 0, "S")
    return 1./gr.Eval(x[0]) # linear interpolation sufficient


def foo(x):
    #return gr.Eval(x[0], 0, "S")
    return 0.01 * x[0] # linear interpolation sufficient

def getGraphFromData(datafile, structure):

    f = open(datafile,'r')
    lines = f.read().splitlines()
    V = []
    meas = []
    err = []
    current = []
    for line in lines:
        if '#' in line: continue
        data = line.split()
        V.append(-1.*float(data[0]))
        if 'MOS' in structure:
            meas.append(1.e12*float(data[-3])) # capacitance in pF (usually order of 1e-10 F in the file)
            err.append(0)
            current.append(1.e06*float(data[-1])) # current in microA
        elif 'GCD' in structure:
            meas.append(-1.e09*float(data[2])) # in nA
            err.append(1e09*float(data[3]))
    tge = ROOT.TGraphErrors()
    for i, v in enumerate(V):
        # do not read above self.GCD_maxV for GCD
        if "GCD" in structure:
            if (v > self.GCD_maxV or meas[i] > 1000.0 or meas[i] < 0.0):
                continue
            #if self.name == "N4789-10_UL" and dose == 1000:
            #    print(f"{self.name}, dose={dose} kGy: meas[{v}] = {meas[i]}")
        tge.SetPoint(i,v,meas[i])
        tge.SetPointError(i,0,err[i])
    yTitle = "MOS capacitance [pF]" if "MOS" in structure else "Diode current [nA]"
    #tge.SetTitle(f"{structure}; Gate voltage [-V]; {yTitle}")
    tge.SetTitle(f"{structure}; Gate voltage [-V]; {yTitle}")
    return tge

def fitVfb(structure, tge, Cox, outdir, step, skipPlot=False):

    plotdir = f"{outdir}/fits/"
    createPlotDirAndCopyPhp(plotdir)

    plotnameTag = f"{plotdir}/fit_{structure}_annealingStep{step}"

    tge = cutGraph(tge)

    c = ROOT.TCanvas()
    c.SetRightMargin(0.06)
    c.SetTickx(1)
    c.SetTicky(1)

    slope = getDerivative(tge)
    slope.GetXaxis().SetTitle("Voltage [-V]")
    slope.GetYaxis().SetTitle("#Delta C / #Delta V [pF/V]")
    slope.GetXaxis().SetTitleSize(0.045)
    slope.GetYaxis().SetTitleSize(0.045)
    slope.GetYaxis().SetTitleOffset(0.99)
    slope.SetMarkerStyle(20)
    slope.Draw("apl")
    if not skipPlot:
        nameNoExt = f"{plotnameTag}_firstDerivative"
        for ext in ["png"]:
            c.SaveAs(f"{nameNoExt}.{ext}")

    minC = min(list(tge.GetY()))
    maxC = max(list(tge.GetY()))
    maxV = max(list(tge.GetX()))
    diff = maxC-minC

    # get x points around which the slope graph has a maximum, which is were it mostly stabilizes
    # this can't work if there are too few points
    maxSlope = max(slope.GetY())
    threshold = 0.8 * maxSlope
    threshold_alt = 0.7 * maxSlope
    low_ramp = 0.0
    high_ramp = 0.0
    low_ramp_alt = 0.0
    high_ramp_alt = 0.0
    xflatSlope = [slope.GetPointX(i) for i in range(slope.GetN()) if slope.GetPointY(i) > threshold]
    if len(xflatSlope) > 3:
        low_ramp = min(xflatSlope)
        high_ramp = max(xflatSlope)
        xflatSlope = [slope.GetPointX(i) for i in range(slope.GetN()) if slope.GetPointY(i) > threshold_alt]
        low_ramp_alt = min(xflatSlope)
        high_ramp_alt = max(xflatSlope)
    else:
        low_ramp  = findX(minC + 0.5 * diff, tge)
        high_ramp = findX(minC + 0.9 * diff, tge)
        # alternative fit with modified range, to estimate uncertainty on Vfb
        low_ramp_alt  = findX(minC + 0.3 * diff, tge)
        high_ramp_alt = findX(minC + 0.7 * diff, tge)

    ramp = ROOT.TF1('ramp', 'pol1(0)', low_ramp, high_ramp)
    tge.Fit(ramp, 'q', '', low_ramp, high_ramp)
    ramp_alt = ROOT.TF1('ramp_alt', 'pol1(0)', low_ramp_alt, high_ramp_alt) # may add + 10 to high_ramp_alt for better visualization, but still fit in range up to high_ramp_alt     
    ramp_alt.SetLineColor(ROOT.kOrange+2)
    tge.Fit(ramp_alt, 'q0+', '', low_ramp_alt, high_ramp_alt)

    Vfb = -1. * (ramp.GetParameter(0) - Cox) / ramp.GetParameter(1)
    Vfb_alt = -1. * (ramp_alt.GetParameter(0) - Cox) / ramp_alt.GetParameter(1)
    # parameters are partially anticorrelated, so raise one and decrease the other (should use covariance matrix, though)
    deltaVfb = 0.5 * abs(Vfb_alt - Vfb)
    # fit errors negligible
    #Vfb_paramVarUp = -1. * (ramp.GetParameter(0)+ramp.GetParError(0) - Cox) / (ramp.GetParameter(1)-ramp.GetParError(1))
    #Vfb_paramVarDown = -1. * (ramp.GetParameter(0)-ramp.GetParError(0) - Cox) / (ramp.GetParameter(1)+ramp.GetParError(1))
    #deltaVfb = 0.5 * (abs(Vfb_paramVarUp - Vfb) + abs(Vfb_paramVarDown - Vfb))
    Vfb_bestEstimate = 0.5 * (Vfb + Vfb_alt)

    tge.SetMarkerStyle(20)
    tge.SetMarkerSize(0.3)
    tge.GetXaxis().SetTitle("Voltage [-V]")
    tge.GetYaxis().SetTitle("MOS capacitance [pF]")
    tge.GetXaxis().SetTitleSize(0.045)
    tge.GetYaxis().SetTitleSize(0.045)
    tge.GetYaxis().SetTitleOffset(0.99)
    tge.Draw('ap')
    ramp.Draw('l same')
    ramp_alt.Draw('l same')

    ramp_ext = ROOT.TF1('ramp_ext', 'pol1(0)', high_ramp, Vfb+5)
    plat_ext = ROOT.TF1('plat_ext', 'pol0(0)', Vfb-10, maxV)
    ramp_ext.SetParameters(ramp.GetParameter(0), ramp.GetParameter(1))
    plat_ext.SetParameter(0, Cox)
    ramp_ext.SetLineColor(ROOT.kBlue)
    plat_ext.SetLineColor(ROOT.kBlue)

    ramp_ext_alt = ROOT.TF1('ramp_ext_alt', 'pol1(0)', high_ramp_alt, Vfb_alt+5)
    ramp_ext_alt.SetParameters(ramp_alt.GetParameter(0), ramp_alt.GetParameter(1))
    ramp_ext_alt.SetLineColor(ROOT.kAzure+1)
    
    ramp_ext.Draw('l same')
    plat_ext.Draw('l same')
    ramp_ext_alt.Draw('l same')

    l = ROOT.TLine(Vfb, minC+0.2*diff, Vfb, maxC)
    l.SetLineColor(ROOT.kGreen+1)
    l.Draw('l same')
    l_alt = ROOT.TLine(Vfb_alt, minC+0.2*diff, Vfb_alt, maxC)
    l_alt.SetLineColor(ROOT.kGreen+2)
    l_alt.Draw('l same')

    lat = ROOT.TLatex()
    lat.SetNDC();
    lat.SetTextFont(42)
    lat.SetTextSize(0.04)
    #lat.DrawLatex(0.45, 0.2, f"dose = {dose} kGy")
    lat.DrawLatex(0.45, 0.15, "flat-band voltage = {:0.1f} +/- {:0.1f} V".format(Vfb_bestEstimate, deltaVfb))

    if not skipPlot:
        nameNoExt = f"{plotnameTag}"
        for ext in ["png", "pdf"]:
            c.SaveAs(f"{nameNoExt}.{ext}")
            
    return (Vfb_bestEstimate, deltaVfb)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",  type=str, 
                        default="/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/", help="Input folder")
    parser.add_argument("-o", "--outdir", type=str, 
                        default="./allplots", help="Output folder for plots")
    parser.add_argument("-c", "--configdir", type=str, 
                        default=None, help="Folder with configuration files")
    parser.add_argument("-t", "--temperature-log", dest="temperatureLog", required=True, type=str, 
                        default="", help="Folder with temperature data (inside args.indir)") # e.g. temperatureLogs/temperatureLog-2022-02-25-17-25-55.dat
    # skip some of last folders if needed
    parser.add_argument("--last-step", dest="lastStep", type=int, 
                        default=-1, help="Last annealing step to include (starts from 0)")
    # this is mainly to get data for irradiation steps
    parser.add_argument("--max-dose", dest="maxDose", type=int, 
                        default=0, help="Maximum dose to use for all samples, negative means to use all (usually only needs dose 0 to get Cox parameter for MOS)")
    # samples
    parser.add_argument("-s", "--sample", type=str, default='N4789-12_UL', help="Sample to be used")
    parser.add_argument("--structures", nargs="+", default=['MOShalf','MOS2000'], help="List of structures to use")
    # others
    parser.add_argument("--skip-points", dest="skipPoints", nargs=2, action='append', type=int, default=[], help="Points to exclude from fit and plots, passed as comma separated pair like 'in,fin'. Can specify multiple times")
    parser.add_argument("-r", "--add-ratio", dest="addRatio", type=str, 
                        default=None, help="Add ratio plots using the sample passed here at denominator")
    parser.add_argument("--skip-plot", dest="skipPlot", action="store_true", default=False, help="Skip plots of C vs V, saves time")
    parser.add_argument("--alt-fit", dest="doAltFit", action="store_true", default=False, help="Do alternate fits splitting points in two ranges (hardcoded inside for now)")

    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    leftMargin = 0.15
    rightMargin = 0.04

    c = ROOT.TCanvas("c", "")
    c.SetFillColor(0)
    c.SetGrid()
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(leftMargin)
    c.SetRightMargin(rightMargin)
    c.SetTickx(1)
    c.SetTicky(1)

    skipPoints = args.skipPoints
    #print(skipPoints)

    #tgeTemperature = ROOT.TGraphErrors()
    #tgeTemperature.SetName(f"temperature")
    firstTime = None
    temperatureFile = f"{args.indir}/{args.temperatureLog}"    
    with open(temperatureFile) as tf:
        for line in tf:
            if "pt1000" in line: 
                continue
            #print(line.split(" "))
            try:
                date,time_temperature = [x.strip() for x in line.split(" ")]
            except:
                #print(line.split(" "))
                continue
            #print(date)
            time,temperature = list(x.strip() for x in time_temperature.split("\t"))[:2]
            #print(time,temperature)
            time = time.rstrip(":")
            year,month,day = map(int, date.split("-"))
            hour,minute,second = map(int, time.split(":"))
            thistime = datetime.datetime(year, month, day, hour, minute, second)
            #print(year,month,day,hour,minute,second,temperature)            
            #print(thistime)
            if firstTime == None:
                firstTime = thistime
                print("-"*30)
                print("CHECK FIRST TIME")
                print(firstTime)
                print("-"*30)
            #diffminutes = (thistime - firstTime).total_seconds() / 60.
            #tgeTemperature.SetPoint(tgeTemperature.GetN(), diffminutes, float(temperature))
            diffhours = (thistime - firstTime).total_seconds() / 3600.
            tgeTemperature.SetPoint(tgeTemperature.GetN(), diffhours, float(temperature))
        
    temperatureValues = list(tgeTemperature.GetY())
    averageTemperature = sum(temperatureValues) / len(temperatureValues)
    averageTemperatureKelvin = averageTemperature + 273.15
    minTemp = min(temperatureValues)
    maxTemp = max(temperatureValues)
    print(" "*30)
    print(f"minTemp = {minTemp}      maxTemp = {maxTemp}      aveTemp = {round(averageTemperature,1)} [C]")
    print(" "*30)
    xmin = tgeTemperature.GetPointX(0)
    xmax = tgeTemperature.GetPointX(tgeTemperature.GetN()-1)
    deltax = xmax - xmin

    sample = args.sample

    # mostly to get Cox, only needs dose = 0
    siSeSa = siliconSensorSample(sample, -20, options=args, noGCD=True)

    outdir = f"{args.outdir}/{sample}/"
    createPlotDirAndCopyPhp(outdir)

    regexp = f"{sample}_.*_[0-9]+kGy_annealingStep[0-9]+"
    regMatch = re.compile(regexp)
    folders = [f for f in os.listdir(args.indir) if regMatch.match(f)]
    #print(folders)
    #steps = sorted(list(filter(lambda x: int(x.split("_annealingStep")[-1]), folders)))
    annealingSteps = {}
    for f in folders:
        annStep = int(str(f.split("_")[-1]).lstrip("annealingStep"))
        subfolder = str(os.listdir(args.indir + "/" + f)[0])
        timestamp = subfolder.split("_")[1:]
        path = f"{args.indir}/{f}/{subfolder}/" 
        annealingSteps[annStep] = {"timestamp": timestamp,
                                   "path"     :  path}

    for structure in args.structures:
        tgeVsTime = ROOT.TGraphErrors()
        tgeVsTime.SetName(f"sample_{structure}_annealing")
        c.cd()
        for k in sorted(annealingSteps.keys()):
            skip = False
            if args.lastStep >= 0 and k > args.lastStep:
                skip = True
            for kskip in skipPoints:
                if k <= kskip[1] and k >= kskip[0]:
                    skip = True
            if skip:
                print(f"Skipping step {k}")
                continue
            #print(f"{k} -> {annealingSteps[k]}")
            datafiles = list(filter(lambda x: x.endswith(".dat"), os.listdir(annealingSteps[k]["path"])))
            Cox = siSeSa.getCox(structure)
            onefile = list(x for x in datafiles if structure in x)[0] 
            datafile = annealingSteps[k]["path"] + "/" + onefile
            #print(datafile)
            graph = getGraphFromData(datafile, structure)
            Vfb, deltaVfb = fitVfb(structure, graph, Cox, outdir, k, skipPlot=args.skipPlot)
            date = annealingSteps[k]["timestamp"][0]
            year = int(date[0:4])
            month = int(date[4:6])
            day = int(date[6:8])
            time = annealingSteps[k]["timestamp"][1]
            hour = int(time[0:2])
            minute = int(time[2:4])
            second = int(time[4:6])
            thistime = datetime.datetime(year, month, day, hour, minute, second)
            #diffminutes = (thistime - firstTime).total_seconds() / 60.
            diffhours = (thistime - firstTime).total_seconds() / 3600.
            #print(f"{k} -> dt = {diffminutes}")
            tgeVsTime.SetPoint(tgeVsTime.GetN(), diffhours, Vfb)
            tgeVsTime.SetPointError(tgeVsTime.GetN()-1, 0, deltaVfb)
            #tgeVsTime.SetPointError(tgeVsTime.GetN()-1, 0, 0.01*Vfb) # dummy in case uncertainties are ill-defined, usually they are around 1%

        graphValues = list(tgeVsTime.GetY())
        minVal = min(graphValues)
        maxVal = max(graphValues)
        valDiff = maxVal - minVal
        scale = 0.5 * valDiff / (maxTemp - minTemp) # scale size of temperature band to 50% of other graph band
        tgeTemperature_scaled = tgeTemperature.Clone(f"tgeTemperature_scaled_{structure}")
        for i in range(tgeTemperature.GetN()):
            yval = tgeTemperature.GetPointY(i)
            tgeTemperature_scaled.SetPointY(i, (yval - minTemp) * scale + minVal + 0.25 * valDiff)

        maxTempScaled = minVal + 0.75 * valDiff
        minTempScaled = minVal + 0.25 * valDiff
        horizLine = ROOT.TLine()
        horizLine.SetLineColor(ROOT.kGray+3)
        horizLine.SetLineWidth(2)
        horizLine.SetLineStyle(2)

        yTitle = "GCD current [nA]" if structure == "GCD" else "Flat-band voltage [-V]"
        tgeVsTime.SetTitle(f"Sample: {sast.getSampleAttribute(sample, 'leg')}   Structure: {structure}; Time [hours]; {yTitle}")
        tgeVsTime.SetLineColor(ROOT.kAzure+2)
        tgeVsTime.SetLineWidth(2)
        tgeVsTime.SetMarkerStyle(7)
        tgeVsTime.Draw('ap')
        tgeVsTime.GetXaxis().SetTitleSize(0.045)
        tgeVsTime.GetXaxis().SetRangeUser(-50, 1.15*tgeVsTime.GetPointX(tgeVsTime.GetN()-1))
        tgeVsTime.GetYaxis().SetTitleSize(0.045)
        tgeVsTime.GetYaxis().SetTitleOffset(0.99)
        nomifit = ROOT.TF1("nomifit", "[0]*((1+x/[1])^([2]))", tgeVsTime.GetPointX(0), tgeVsTime.GetPointX(tgeVsTime.GetN()-1))
        #nomifit = ROOT.TF1("nomifit", "[0]*((1-[1]*exp([3]/%s)*x)^([2]))" % round(averageTemperatureKelvin,2), tgeVsTime.GetPointX(0), tgeVsTime.GetPointX(tgeVsTime.GetN()-1))
        if structure == "MOS2000":
            nomifit.SetParameters(tgeVsTime.GetPointY(0), 1.0, -0.07)
        else:
            nomifit.SetParameters(tgeVsTime.GetPointY(0), 1.0, -0.1)
            #nomifit.SetParameters(102, -1.0, -0.1, 5000.0)
        nomifit.SetLineColor(ROOT.kOrange+2)
        tgeVsTime.Fit("nomifit","SEMR+ EX0") 

        # 
        altfit = []
        nAltFitPointEdges = [0, 21, tgeVsTime.GetN()-1]
        altfitColors = [ROOT.kRed+1, ROOT.kGreen+2, ROOT.kCyan+1, ROOT.kMagenta, ROOT.kViolet, ROOT.kGray+2]
        if args.doAltFit:
            for ialt in range(len(nAltFitPointEdges) - 1):
                ipLow = nAltFitPointEdges[ialt]
                ipHigh = nAltFitPointEdges[ialt+1]
                altfit.append(ROOT.TF1(f"altfit_{ialt}", "[0]*((1+x/[1])^([2]))", tgeVsTime.GetPointX(ipLow), tgeVsTime.GetPointX(ipHigh)))
                # set first parameter to value at t = 0 because the fit still uses "t" and not "t - tgeVsTime.GetPointY(ipLow)"
                altfit[ialt].SetParameters(tgeVsTime.GetPointY(0), nomifit.GetParameter(1), nomifit.GetParameter(2))
                altfit[ialt].SetLineColor(ROOT.kBlack if ialt >= len(altfitColors) else altfitColors[ialt])
                print("-"*30)
                print(" "*30)
                print(f"Fit between points ix = {ipLow} and ix = {ipHigh}")
                tgeVsTime.Fit(f"altfit_{ialt}","SEMR+ EX0") 
                print("-"*30)


        # evaluate asymptoticvalue of time for which Vfb < 3 V (approximately the original value)
        Vfb_asym = 3.0
        Vextrap = tgeVsTime.GetPointY(tgeVsTime.GetN()-1)
        timeExtrap = tgeVsTime.GetPointX(tgeVsTime.GetN()-1)
        timeStep = 24.0 # in hours, so 1 day
        i = 0
        while Vextrap > Vfb_asym:
            timeExtrap += timeStep
            Vextrap = nomifit.Eval(timeExtrap)
            i += 1
            if i > 365: break # after one year stop this loop 
        print("-"*30)
        print(f"{structure}: Vextrap = {Vextrap} V @ time = {round(timeExtrap,1)} hours ({round(timeExtrap/24.0,1)} days or {round(timeExtrap/168.0,1)} weeks) ")
        print("-"*30)

        tgeTemperature_scaled.SetLineColor(ROOT.kYellow+2)
        tgeTemperature_scaled.SetLineWidth(2)
        tgeTemperature_scaled.SetLineStyle(3)
        tgeTemperature_scaled.Draw("l same")

        horizLine.DrawLine(xmin, maxTempScaled, xmax, maxTempScaled)
        horizLine.DrawLine(xmin, minTempScaled, xmax, minTempScaled)
        lat = ROOT.TLatex()
        #lat.SetNDC();
        lat.SetTextFont(42)
        lat.SetTextSize(0.04)
        lat.DrawLatex(xmax - 0.15 * deltax, maxTempScaled - 0.06 * valDiff, "T_{max} = %.1f C" % maxTemp)
        lat.DrawLatex(xmax - 0.15 * deltax, minTempScaled - 0.06 * valDiff, "T_{min} = %.1f C" % minTemp)

        leg = ROOT.TLegend(0.35, 0.7, 0.96, 0.9)
        leg.SetNColumns(2)
        leg.SetFillColor(0)
        #leg.SetFillStyle(0)
        #leg.SetBorderSize(0)
        leg.AddEntry(tgeVsTime, "Measurement", "PLE")
        leg.AddEntry(tgeTemperature_scaled, "Temperature", "L")
        #leg.AddEntry(expfit, "Fit: a+b#upoint e^{-t/c}", "L")
        leg.AddEntry(nomifit, "Fit: a #upoint(1 #plus t / b)^{c}", "L")
        leg.Draw("same")

        cName = f"{outdir}/summary_annealingVsTime_{sample}_{structure}"
        for ext in ["png", "pdf"]:
            c.SaveAs(f"{cName}.{ext}")

        # compare temperature with ratio of graph/fit
        # for ratio and comparison to temperature
        canvas = ROOT.TCanvas("canvas","",800,800)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetGridx(1)
        canvas.SetGridy(1)
        canvas.SetLeftMargin(leftMargin)
        canvas.SetRightMargin(rightMargin)
        canvas.cd()
        canvas.SetBottomMargin(0.5)
        pad2 = ROOT.TPad("pad2","pad2",0.0,0.0,1.0,0.95)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.SetTopMargin(0.5)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetRightMargin(rightMargin)
        pad2.SetFillColor(0)
        pad2.SetGridx(1)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)
        ratio = ROOT.TGraphErrors()
        ratio.SetName(f"measurementOverFit_{structure}")
        ratio.SetTitle(f"Sample: {sast.getSampleAttribute(sample, 'leg')}   Structure: {structure}")
        for i in range(tgeVsTime.GetN()):
            xval = tgeVsTime.GetPointX(i)
            fitval = nomifit.Eval(xval)
            ratio.SetPoint(i, xval, tgeVsTime.GetPointY(i)/fitval)
            ratio.SetPointError(i, 0, tgeVsTime.GetErrorY(i)/fitval)
        canvas.cd()
        #rangeXaxis = 1.05 * (max(tgeVsTime.GetPointX(tgeVsTime.GetN()-1), xmax) - tgeTemperature.GetPointX(0))
        #hframe = ROOT.TH1D("hframe", "", 1, tgeTemperature.GetPointX(0), rangeXaxis)
        #hframe.SetMarkerSize(0)
        ratio.SetMarkerStyle(20)
        ratio.GetYaxis().SetTitle("Measured V_{fb} / fit")
        ratio.GetXaxis().SetLabelSize(0)
        ratio.GetXaxis().SetTitle("")
        #hframe.GetXaxis().SetRangeUser(tgeTemperature.GetPointX(0), rangeXaxis)
        #minyframe = min(list(ratio.GetY()))
        #maxyframe = max(list(ratio.GetY()))
        #diff = maxyframe - minyframe 
        #hframe.GetYaxis().SetRangeUser(minyframe - 0.05*diff, maxyframe + 0.05*diff)
        #hframe.Draw("HIST")
        #ratio.Draw("apl SAME")
        ratio.Draw("apl")
        pad2.Draw()
        pad2.cd()
        tgeTemperature.GetXaxis().SetTitle("Time [hours]")
        tgeTemperature.GetYaxis().SetTitle("Temperature [C]")
        tgeTemperature.GetYaxis().SetTitleSize(0.05)
        tgeTemperature.GetYaxis().SetLabelSize(0.04)
        #tgeTemperature.GetXaxis().SetRangeUser(ratio.GetPointX(0), rangeXaxis)
        #deltaT = maxTemp - minTemp 
        #hframe.GetYaxis().SetRangeUser(minTemp - 0.05*deltaT, maxTemp + 0.05*deltaT)
        #hframe.Draw("HIST")
        #tgeTemperature.Draw("apl SAME")
        tgeTemperature.Draw("apl")

        cName = f"{outdir}/measurementOverFit_{sample}_{structure}"
        for ext in ["png", "pdf"]:
            canvas.SaveAs(f"{cName}.{ext}")
            
        # save graphs in root file, including fit function, why not
        of = safeOpenFile(f"{outdir}/{sample}_{structure}.root", mode="RECREATE")
        of.cd()
        tgeTemperature.Write()
        ratio.Write()
        nomifit.Write()
        tgeVsTime.Write()
        tgeTemperature_scaled.Write()
        of.Close()
