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

class siliconSensorSample:
    def __init__(self, name, temperature, options=None, title=None): 
        self.options = options
        self.name = name
        self.typeName = self.getTypeName()
        self.typeNameGoodString = self.typeName.replace(' ','_').replace('#','v') 
        self.title = name if title == None else title
        self.datapath = self.options.indir # usually /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/logs_obelix_setup/
        self.temperature = temperature # temperature for measurement, usually -20 C
        self.structures = self.options.structures # usuallly ['MOShalf','MOS2000', 'GCD']
        self.configdir = self.options.configdir
        self.plotdir = f"{self.options.outdir}/{self.typeNameGoodString}_{self.name}"
        self.GCD_cuts = []
        self.GCD_maxV = 85.0
        if self.configdir:
            self.readGCDcuts()
        createPlotDirAndCopyPhp(self.plotdir)

        self.doses = {s: self.getAllDoses(s) for s in self.structures}  # are the doses always the same for all structures?
        #self.graphs = dict.fromkeys(self.structures, None) # will have a graph per structure
        self.graphsVsDose = {s : None for s in self.structures} # will have a graph per structure, with Vfb vs dose
        self.graphsSingleDose = {s : {d: None for d in self.doses[s]} for s in self.structures}
        self.MOScurrentSingleDose = {}
        if "EPI" in self.typeName:
            self.MOScurrentSingleDose = {s : {d: None for d in self.doses[s]} for s in self.structures if "MOS" in s}
        self.Cox = {s : None for s in self.structures if "MOS" in s} # to be filled later once graphs exist

        self.readData()
        self.plotData()
        if "EPI" in self.typeName:
            self.plotMOScurrent()

    def getName(self):
        return self.name

    def getTypeName(self):
        return getSampleTypeFromName(self.name)

    def getGraphVsDose(self, structure):
        return self.graphsVsDose[structure]

    def getGraphSingleDose(self, structure, dose):
        return self.graphsSingleDose[structure][dose]

# should this be a method of siliconSensorSample?
    def getPath(self, mainpath, sample, structure, dose, temperatureFloat):
        temperature = str(temperatureFloat).replace('-','m')
        prefix = "cv" if "MOS" in structure else "iv" if "GCD" in structure else "unknown" 
        path = f"{mainpath}/{sample}_{temperature}C_{dose}kGy"
        path = glob(f"{path}/*/")[-1]
        path += f"{prefix}_{sample}_{temperature}C_{dose}kGy_{structure}.dat"
        if not os.path.exists(path):
            print(f"Error: I did not find this file {path}")
            quit()
        return path

    def getAllDoses(self, structure):
        temperature = str(self.temperature).replace('-','m')
        prefix = "cv" if "MOS" in structure else "iv" if "GCD" in structure else "unknown" 
        regexp = f"{self.name}_{temperature}C_[0-9]+kGy"
        #print(regexp)
        regMatch = re.compile(regexp)
        folders = [f for f in os.listdir(self.datapath) if regMatch.match(f)]
        #print(folders)
        doses = []
        for f in folders:
            #print(f"{structure}: {f}")
            dose = str(f.split("_")[-1]).rstrip("kGy")
            doses.append(int(dose))
        doses = sorted(doses)
        #print(f"{structure}: {doses}")
        return doses

    def readGCDcuts(self):
        if not self.configdir:
            return
        f = open(f'{self.configdir}/cut_GCD_curve.txt','r')
        lines = f.read().splitlines()
        f.close()
        lines = list(filter(lambda a: not a.startswith('#') and not a.replace(' ','') == '', lines))
        self.GCD_cuts = [l.replace(' ','').split('&') for l in lines]
        for l in self.GCD_cuts:
            if len(l) != 4 and len(l) != 6:
                print(f"ERROR in file cut_GCD_curve.txt: wrong number of arguments in entry below\n {l}")
                sys.exit()
        return

    def getGCDcuts(self, dose):
        low = high = gmin = gmax = -999
        for l in self.GCD_cuts:
            if self.name != l[0]: continue
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

    def cutGCDcurve(self, tge, cut):
        if list(tge.GetX())[-1] > cut:
            tge.RemovePoint(tge.GetN()-1)
        if list(tge.GetX())[-1] > cut:
            tge = self.cutGCDcurve(tge,cut)
        tge.SetMaximum(max(list(tge.GetY()))*1.1)
        return tge

    def cutGCDcurveLow(self, tge, cut):
        if list(tge.GetX())[0] < cut:
            tge.RemovePoint(0)
        if list(tge.GetX())[0] < cut:
            tge = self.cutGCDcurveLow(tge,cut)
        tge.SetMaximum(max(list(tge.GetY()))*1.1)
        return tge

    def removeBadPoints(self, tge, threshold):
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
            tge = self.removeBadPoints(tge, threshold)
        return tge

    def findApproxDepletion(self, dose, tge):
        X = list(tge.GetX())
        Y = list(tge.GetY())
        eY = list(tge.GetEY())
        baseline = min(Y)
        Y = [y-baseline for y in Y]
        for i,y in enumerate(Y):
            if y > max(Y)*0.3: 
                break

        deltaX = 10 if dose < 40 else 20
        xL = X[i] - deltaX
        xH = X[i] + deltaX
        # for now try without next line
        #xL, xH = getGCDrange(self.name, dose, xL, xH)
        if self.name == "N4790-1_LR" and dose >= 10:
            xL = 20.0
            xH = 50.0 if dose <= 40 else 60.0
        if self.name == "N4790-13_LR" and dose >= 10 and dose < 70:
            xL = 16.0
            xH = 45.0

        xH = min(85, xH) # better not to go above about self.GCD_maxV V

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
        return [xm, ym+baseline, xM, yM+baseline, (em**2+eM**2)**0.5]

    def readData(self):

        for structure in self.structures:
            for dose in self.doses[structure]:
                datafile = self.getPath(self.datapath, self.name, structure, dose, self.temperature)
                #print(datafile)
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
                        meas.append(1.e09*float(data[-3]))
                        err.append(0)
                        current.append(1.e06*float(data[-1])) # current in microA
                    elif 'GCD' in structure:
                        meas.append(-1.e09*float(data[2]))
                        err.append(1e09*float(data[3]))
                        
                tge = ROOT.TGraphErrors()
                if len(self.MOScurrentSingleDose.keys()) and "MOS" in structure:
                    self.MOScurrentSingleDose[structure][dose] = ROOT.TGraphErrors()
                    self.MOScurrentSingleDose[structure][dose].SetTitle(f"{self.name} {structure}; Gate voltage [-V]; MOS current [-#muA]")
                for i, v in enumerate(V):
                    # do not read above self.GCD_maxV for GCD
                    if "GCD" in structure:
                        if (v > self.GCD_maxV or meas[i] > 1000.0 or meas[i] < 0.0):
                            continue
                        #if self.name == "N4789-10_UL" and dose == 1000:
                        #    print(f"{self.name}, dose={dose} kGy: meas[{v}] = {meas[i]}")
                    tge.SetPoint(i,v,meas[i])
                    tge.SetPointError(i,0,err[i])
                if len(self.MOScurrentSingleDose.keys()) and "MOS" in structure:
                    for i, v in enumerate(V):
                        self.MOScurrentSingleDose[structure][dose].SetPoint(i,v,-1.0*current[i]) # current is negative, let's make it positive
                yTitle = "MOS capacitance [pF]" if "MOS" in structure else "Diode current [nA]"
                #tge.SetTitle(f"{structure}; Gate voltage [-V]; {yTitle}")
                self.graphsSingleDose[structure][dose] = tge
                self.graphsSingleDose[structure][dose].SetTitle(f"{structure}; Gate voltage [-V]; {yTitle}")

            if "MOS" in structure:
                # store Cox for later usage, no need to recompute it everytime
                self.getCox(structure)

        return

    def getCox(self, structure):
        # use 0 kGy dose to get the capacitance plateau
        dose = 0
        tge = self.graphsSingleDose[structure][dose]
        if tge is None:
            print(f"Warning: graph is None for {structure} and dose = {dose}")
            quit()
        Cox = max(list(tge.GetY()))
        #print(f"Cox[{structure}] = {Cox}")
        self.Cox[structure] = Cox
        return Cox

    # to get surface velocity from the GCD current
    # def calculate_GCD_parameters(self, current):
    #     current *= 1E-09 # in A
    #     A = cnst.A_GCD
    #     # if self.name == '3009_LR' or self.name == '3010_LR':
    #     #     A = cnst.A_GCD_6inch
    #     # elif '_SE_' in self.name:
    #     #     A = cnst.A_GCD_tracker
    #     A *= 1E-6 # in m^2
    #     ni = cnst.ni * 1E+06 # in m^-3
    #     J = current / (cnst.q*ni*A) # in m/s
    #     J *= 100 # in cm/s
    #     return J 

    def plotData(self):
        rfName = f"{self.plotdir}/summaryVsDose_{self.name}"
        tf = safeOpenFile(f"{rfName}.root", mode='recreate')
        for structure in self.structures:
            tgeVsDose = ROOT.TGraphErrors()
            tgeVsDose.SetName(f"{self.name}_{structure}")
            for dose in self.doses[structure]:
                # check if maximum voltage was already reached at this dose, in that case the graph is not reliable
                # for MOS will evaluate maxY wrt Cox, which should be the expected maximum capacitance to be reached
                # if the maxY is below Cox by more than 5% of the difference between Cox and minY 
                ylist = list(self.graphsSingleDose[structure][dose].GetY())
                maxY = max(ylist)
                minY = min(ylist)
                if "MOS" in structure:
                    diff = self.Cox[structure] - minY
                    Vfb, deltaVfb = self.fitVfb(structure, dose)
                    # set point in graph vs dose only if the ramp up arrived close to the top 
                    if (self.Cox[structure] - maxY) < (0.15 * diff):
                        tgeVsDose.SetPoint(tgeVsDose.GetN(), dose, Vfb)
                        #tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, deltaVfb)
                        tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, 0) # no error on Y for now
                else:
                    [current, err] = self.getGCDcurrent(dose)
                    # do not use the point for dose = 0 for GCD
                    if dose:
                        tgeVsDose.SetPoint(tgeVsDose.GetN(), dose, current)
                        tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, err)
                    # for surface velocity, but not needed for now
                    #J = self.calculate_GCD_parameters(current)                    
                    #tgeVsDose.SetPoint(tgeVsDose.GetN(), dose, J)
                    #tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, err*J/current)

            self.graphsVsDose[structure] = tgeVsDose

            c = ROOT.TCanvas()
            c.SetRightMargin(0.06)
            c.SetTickx(1)
            c.SetTicky(1)
            yTitle = "GCD current [nA]" if structure == "GCD" else "Flat-band voltage [-V]"
            self.graphsVsDose[structure].SetTitle(f"{self.name} {structure}; Dose [kGy]; {yTitle}")
            self.graphsVsDose[structure].SetMarkerStyle(7)
            self.graphsVsDose[structure].Draw('apl')
            self.graphsVsDose[structure].GetXaxis().SetTitleSize(0.045)
            self.graphsVsDose[structure].GetXaxis().SetRangeUser(-20, 1.15*tgeVsDose.GetPointX(tgeVsDose.GetN()-1))
            self.graphsVsDose[structure].GetYaxis().SetTitleSize(0.045)
            self.graphsVsDose[structure].GetYaxis().SetTitleOffset(0.99)
            cName = f"{self.plotdir}/summaryVsDose_{self.name}_{structure}"
            for ext in ["png", "pdf"]:
                c.SaveAs(f"{cName}.{ext}")
            self.graphsVsDose[structure].Write(f"{self.name}_{structure}")
            c.Clear()

        tf.Close()

        return

    def plotMOScurrent(self):

        c = ROOT.TCanvas()
        c.SetTickx(1)
        c.SetTicky(1)
        c.cd()
        c.SetFillColor(0)
        c.SetGrid()
        c.SetLeftMargin(0.14)
        c.SetRightMargin(0.06)
        c.cd()

        ROOT.gStyle.SetPalette(87)
        for structure in self.MOScurrentSingleDose.keys():
            minX = 0
            maxX = 0
            maxY = 0
            minY = 0
            sortedDoses = sorted(list(self.MOScurrentSingleDose[structure].keys()))
            graphs = {d : self.MOScurrentSingleDose[structure][d] for d in sortedDoses} 
            graphsAverage = {d : ROOT.TGraphErrors() for d in sortedDoses} 
            for dose in graphs.keys():
                gr = graphs[dose]
                # now make averages 
                nPtsPerAverage = 10 if dose < 10 else 20
                allXvalues = sorted(list(gr.GetX()))
                sumY = 0.0
                sumX = 0.0
                nPts = 0
                for i,x in enumerate(allXvalues):
                    sumY += gr.Eval(x)
                    sumX += x
                    nPts += 1
                    if nPts == nPtsPerAverage:
                        graphsAverage[dose].SetPoint(graphsAverage[dose].GetN(), sumX/nPts, sumY/nPts)
                        sumY = 0.0
                        sumX = 0.0
                        nPts = 0
                if nPts > 0:
                    graphsAverage[dose].SetPoint(graphsAverage[dose].GetN(), sumX/nPts, sumY/nPts)
                #allYvalues = list(gr.GetY())
                allYvalues = list(graphsAverage[dose].GetY())
                maxX = max(maxX, graphsAverage[dose].GetPointX(graphsAverage[dose].GetN()-1))
                maxY = max(maxY, max(allYvalues))
                minY = min(minY, min(allYvalues))
            maxX *= 1.05
            maxY *= 1.8
            minY *= 0.5
            frame = ROOT.TH1D(f"frame_MOScurrent_{self.typeName}_{structure}", f"Sample: {self.typeName}   {structure}", 1, minX, maxX)
            frame.SetMarkerSize(0)
            frame.SetMarkerColor(0)
            frame.SetLineColor(0)
            frame.Draw()
            frame.GetXaxis().SetTitleSize(0.05)
            frame.GetXaxis().SetLabelSize(0.04)
            frame.GetXaxis().SetTitleOffset(0.9)
            frame.GetYaxis().SetTitleOffset(1.3)
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetYaxis().SetLabelSize(0.04)
            frame.GetXaxis().SetTitle("Gate voltage [-V]")
            frame.GetYaxis().SetTitle("MOS current [-#muA]")
            frame.GetXaxis().SetRangeUser(minX,maxX)
            frame.GetYaxis().SetRangeUser(minY,maxY)
            frame.SetStats(0)
            
            leg = ROOT.TLegend(0.2,0.68,0.9,0.88)
            leg.SetFillColorAlpha(0,0.6)
            leg.SetBorderSize(0)
            leg.SetNColumns(4)

            for dose in graphs.keys():
                graphsAverage[dose].Draw("PL PLC PFC SAME")
                leg.AddEntry(graphsAverage[dose], f"{dose} kGy", "PLF")
                #graphs[dose].Draw("PL PLC PFC SAME")
                #leg.AddEntry(graphs[dose], f"{dose} kGy", "PLF")
                #gr.SetLineWidth(2)
            c.RedrawAxis("sameaxis")
            leg.Draw("SAME")

            cName = f"{self.plotdir}/MOScurrentVsDose_{self.name}_{structure}"
            for ext in ["png", "pdf"]:
                c.SaveAs(f"{cName}.{ext}")

        c.Clear()
        return

    def getGCDcurrent(self, dose):
        tge = self.graphsSingleDose["GCD"][dose]
        if tge == None: 
            return

        low, high, gmin, gmax = self.getGCDcuts(dose)
        high = min(self.GCD_maxV, high) # better not to go above 85 V
        for ix,xval in enumerate(list(tge.GetX())):
            if xval > high:
                tge.RemovePoint(ix)

        tge_orig = tge.Clone()    


        if low != -999:
            tge = self.cutGCDcurveLow(tge, low)
        if high != -999:
            tge = self.cutGCDcurve(tge, high)
        if gmin != -999:
            tge.SetMinimum(gmin)
        if gmax != -999:
            tge.SetMaximum(gmax)

        threshold = 0.05
        tge = self.removeBadPoints(tge, threshold)

        xm, ym, xM, yM, err = self.findApproxDepletion(dose, tge)
        deltaI = (max(0.001, yM - ym))

        at = ROOT.TGraph()
        at.SetPoint(0, xm, ym)
        at.SetPoint(1, xM, yM)
        at.SetMarkerStyle(8)
        at.SetMarkerColor(ROOT.kGreen+2)

        c = ROOT.TCanvas()
        tge_orig.SetLineColor(ROOT.kRed+2)
        tge_orig.SetMarkerColor(ROOT.kRed+2)
        
        tge_orig.SetMarkerStyle(20)
        tge_orig.SetMarkerSize(0.3)
        tge_orig.GetXaxis().SetTitle("Voltage [V]")
        tge_orig.GetYaxis().SetTitle("Current [nA]")
        tge_orig.GetXaxis().SetTitleSize(0.045)
        tge_orig.GetYaxis().SetTitleSize(0.045)
        tge_orig.GetYaxis().SetTitleOffset(0.99)

        leg = ROOT.TLegend(0.1, 0.5, 0.5, 0.65)
        leg.SetFillColorAlpha(0,0.6)
        leg.SetBorderSize(0)
        leg.SetNColumns(1)
        
        tge_orig.Draw('apl')
        tge.Draw('pl same')
        at.Draw('p same')
        leg.AddEntry(tge_orig,"Original", "PL")
        leg.AddEntry(tge,     "Polished", "PL")
        leg.AddEntry(at,      "#DeltaI = {:0.2f} #pm {:0.2f}".format(deltaI, err), "P")
        leg.Draw("same")

        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)
        lat.SetTextSize(0.04)
        lat.DrawLatex(0.2, 0.85, f"dose = {dose} kGy")

        for ext in ["png", "pdf"]:
            c.SaveAs(f"{self.plotdir}/fit_{self.name}_GCD_{dose}kGy.{ext}")

        return [deltaI, err]

    
    def fitVfb(self, structure, dose):

        tge = self.graphsSingleDose[structure][dose]
        tge = cutGraph(tge)

        c = ROOT.TCanvas()
        c.SetRightMargin(0.06)
        c.SetTickx(1)
        c.SetTicky(1)

        minC = min(list(tge.GetY()))
        maxC = max(list(tge.GetY()))
        maxV = max(list(tge.GetX()))
        diff = maxC-minC

        low_ramp  = findX(minC + 0.5 * diff, tge)
        high_ramp = findX(minC + 0.9 * diff, tge)
        # alternative fit with modified range, to estimate uncertainty on Vfb
        low_ramp_alt  = findX(minC + 0.3 * diff, tge)
        high_ramp_alt = findX(minC + 0.7 * diff, tge)
        if dose == 0:
            low_ramp  = findX(minC + 0.3 * diff, tge)
            high_ramp = findX(minC + 0.7 * diff, tge)
            low_ramp_alt  = findX(minC + 0.2 * diff, tge)
            high_ramp_alt = findX(minC + 0.6 * diff, tge)

        ramp = ROOT.TF1('ramp', 'pol1(0)', low_ramp, high_ramp)
        tge.Fit(ramp, 'q', '', low_ramp, high_ramp)
        ramp_alt = ROOT.TF1('ramp_alt', 'pol1(0)', low_ramp_alt, high_ramp_alt) # may add + 10 to high_ramp_alt for better visualization, but still fit in range up to high_ramp_alt     
        ramp_alt.SetLineColor(ROOT.kOrange+2)
        tge.Fit(ramp_alt, 'q0+', '', low_ramp_alt, high_ramp_alt)
     
        Vfb = -1. * (ramp.GetParameter(0) - self.Cox[structure]) / ramp.GetParameter(1)
        Vfb_alt = -1. * (ramp_alt.GetParameter(0) - self.Cox[structure]) / ramp_alt.GetParameter(1)
        deltaVfb = abs(Vfb_alt - Vfb)

        tge.SetMarkerStyle(20)
        tge.SetMarkerSize(0.3)
        tge.GetXaxis().SetTitle("Voltage [V]")
        tge.GetYaxis().SetTitle("Capacitance [pF]")
        tge.GetXaxis().SetTitleSize(0.045)
        tge.GetYaxis().SetTitleSize(0.045)
        tge.GetYaxis().SetTitleOffset(0.99)
        tge.Draw('ap')
        ramp.Draw('l same')
        ramp_alt.Draw('l same')

        ramp_ext = ROOT.TF1('ramp_ext', 'pol1(0)', high_ramp, Vfb+5)
        plat_ext = ROOT.TF1('plat_ext', 'pol0(0)', Vfb-10, maxV)
        ramp_ext.SetParameters(ramp.GetParameter(0), ramp.GetParameter(1))
        plat_ext.SetParameter(0, self.Cox[structure])
        ramp_ext.SetLineColor(ROOT.kBlue)
        plat_ext.SetLineColor(ROOT.kBlue)

        ramp_ext.Draw('l same')
        plat_ext.Draw('l same')

        l = ROOT.TLine(Vfb, minC, Vfb, maxC)
        l.SetLineColor(ROOT.kGreen+1)
        l.Draw('l same')

        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)
        lat.SetTextSize(0.04)
        lat.DrawLatex(0.45, 0.2, f"dose = {dose} kGy")
        lat.DrawLatex(0.45, 0.15, "flat-band voltage = {:0.1f} +/- {:0.1f} V".format(Vfb, deltaVfb))

        nameNoExt = f"{self.plotdir}/fit_{self.name}_{structure}_{dose}kGy"
        for ext in ["png", "pdf"]:
            c.SaveAs(f"{nameNoExt}.{ext}")
        return (Vfb, deltaVfb)


# whatever is below this line must be rewritten most likely

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



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",  type=str, 
                        default="/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/logs_obelix_setup/", help="Input folder")
    parser.add_argument("-o", "--outdir", type=str, 
                        default="./allplots", help="Output folder for plots")
    parser.add_argument("-c", "--configdir", type=str, 
                        default=None, help="Folder with configuration files")
    # samples
    parser.add_argument("--samples", nargs="+", default=['N4791-1_LR','N4790-1_UL','N4791-6_UL','N4790-13_LR','N4789-10_UL','N4790-1_LR','N4790-13_UL','N4791-6_LR','N4788-9_LR'], help="List of samples to be used")
    #parser.add_argument("--samples", nargs="+", default=['N4791-1_LR','N4790-1_UL','N4791-6_UL'], help="List of samples to be used")
    parser.add_argument("--structures", nargs="+", default=['MOShalf','MOS2000','GCD'], help="List of structures to use")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    print()
    print(f"Analysing {len(args.samples)} samples")
    print(f"Analysing these structures: {args.structures}")
    print()

    sampleDict = {}
    for sample in args.samples:
        #sampleDict[sample] = siliconSensorSample("N4790-1_UL", args.indir, -20 , plotdir=args.outdir, goodStructures = ['MOShalf','MOS2000'])
        sampleDict[sample] = siliconSensorSample(sample, -20, options=args) #['MOShalf','MOS2000'])

    colors = [ROOT.kBlack, ROOT.kRed+2, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kPink+2, ROOT.kAzure+2, ROOT.kGray+2, ROOT.kSpring+9]
    markers = [ROOT.kFullCircle, ROOT.kOpenCircle, ROOT.kFullTriangleDown, ROOT.kOpenTriangleDown , ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kOpenSquare, ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenSquareDiagonal]
    
    outdirSummary = f"{args.outdir}/summary/"
    createPlotDirAndCopyPhp(outdirSummary)
    
    canvas = ROOT.TCanvas("canvas","",1200,800)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetFillColor(0)
    canvas.SetGrid()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.06)
    canvas.cd()

    #typeNames = [sampleDict[sample].getTypeName() for sample in args.samples]
    #typeNames = sorted(typeNames)

    for structure in args.structures:
        graphs = [sampleDict[sample].getGraphVsDose(structure) for sample in args.samples]
        legendEntries = [sampleDict[sample].getTypeName() for sample in args.samples]
        maxX = 0
        maxY = 0
        minY = 0
        for gr in graphs:
            maxX = max(maxX, gr.GetPointX(gr.GetN()-1))
            maxY = max(maxY, gr.GetPointY(gr.GetN()-1))
            minY = min(minY, gr.GetPointY(gr.GetN()-1))
        maxY *= 1.1
        minX = -20
        useLogX = False
        if maxX > 300:
            useLogX = True
            minX = 0.5 if "GCD" in structure else 0.08
            maxX *= 5.0
            for igr,gr in enumerate(graphs):
                #print(f"{legendEntries[igr]}: setting minX from {gr.GetPointX(0)} to 0.1")
                if "GCD" not in structure:
                    # for GCD we already removed the point at dose=0
                    gr.SetPointX(0, 0.1)

        else:
            maxX *= 1.1

        canvas.SetTitle(structure)
        if "GCD" in structure:
            yTitle = "GCD current [nA]"
            minyLeg = 0.3
            maxyLeg = minyLeg + 0.06 * len(graphs)
            legCoords = f"0.2,{minyLeg},0.5,{maxyLeg}"
            useLogY = False
            #legCoords = f"0.65,{minyLeg},0.95,{maxyLeg}"
        else:
            yTitle = "Flat-band voltage [-V]"
            minyLeg = 0.12
            maxyLeg = minyLeg + 0.06 * len(graphs)
            legCoords = f"0.6,{minyLeg},0.9,{maxyLeg}"
            useLogY = False
            minY = 0
        drawGraphs(graphs, f"Dose [kGy]::{minX},{maxX}", f"{yTitle}::{minY},{maxY}", f"summaryVsDose_compareSamples_{structure}", outdirSummary, 
                   legendEntries=legendEntries, legendCoords=legCoords,
                   vecColors=colors, vecMarkers=markers, passCanvas=canvas, moreText=f"Structure: {structure}", useLogX=useLogX, useLogY=useLogY)
    
        if "GCD" in structure:
            # repeat with log scale on Y axis too
            useLogY = True
            minyLeg = 0.12
            maxyLeg = minyLeg + 0.06 * len(graphs)
            legCoords = f"0.70,{minyLeg},0.95,{maxyLeg}"
            maxY *= 3.0
            minY = max(0.01,minY)
            drawGraphs(graphs, f"Dose [kGy]::{minX},{maxX}", f"{yTitle}::{minY},{maxY}", f"summaryVsDose_compareSamples_{structure}_logY", outdirSummary, 
                       legendEntries=legendEntries, legendCoords=legCoords,
                       vecColors=colors, vecMarkers=markers, passCanvas=canvas, moreText=f"Structure: {structure}", useLogX=useLogX, useLogY=useLogY)


