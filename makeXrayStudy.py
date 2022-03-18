# example

# python3 makeXrayStudy.py -i /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/ -c configs_Nov2021/ -o plots/forMatteoAtMoriond/FZandEPI_noRepetition_TESTRATIO/ --samples 'FZandEPInor' -r 'N4789-10_UL'
# python3 makeXrayStudy.py -i /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/ -c configs_Nov2021/ -o plots/forMatteoAtMoriond/EPItypeC_andAltDose/ --samples 'CandDose' -r 'N4789-10_UL'
# python3 makeXrayStudy.py -i /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/ -c configs_Nov2021/ -o plots/forMatteoAtMoriond/allTypeC_noAltDose_testColor/ --samples 'allCnoDose' -r 'N4789-10_UL'
# python3 makeXrayStudy.py -i /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/ -c configs_Nov2021/ -o plots/forMatteoAtMoriond/onlyFZ_noRepetition/ --samples 'FZnor' -r 'N4791-6_UL'

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
import sampleStyle as sast

from functionsAnnealing import *

import argparse


def getSurfaceVelocity(current):
    I = current * 1E-09 # in A from nA
    A = cnst.A_GCD
    A *= 1E-6 # in m^2
    ni = cnst.ni * 1E+06 # in m^-3
    J = I / (cnst.q*ni*A) # in m/s
    J *= 100 # in cm/s
    return J


def getOxideChargeDensity(V_fb, structure, C_ox):

    if structure not in ['MOShalf', 'MOS2000']:
        print(f"Error in getOxideChargeDensity, unexpected structure {structure}")
    A = cnst.A_MOShalf if structure == 'MOShalf' else cnst.A_MOS2000
    A *= 1E-6 # in m^2
    C_ox -= cnst.approx_openC
    C_ox *= 1E-12 # in F from pF
    phi_s = cnst.Chi + cnst.Eg/2. + cnst.KbTn20*ROOT.log(cnst.NA/cnst.ni)
    phi_ms = cnst.phi_m - phi_s
    N_ox = C_ox / (A*cnst.q) * (phi_ms + V_fb) # 1/m^2
    N_ox *= 1E-04 # in 1/cm^2
    return N_ox


class siliconSensorSample:
    def __init__(self, name, temperature=None, options=None, title=None, noGCD=False): 
        self.options = options
        self.name = name
        self.typeName = self.getTypeName()
        self.typeNameGoodString = self.typeName.replace(' ','_').replace('#','v').replace('/','_').replace('(','_').replace(')','_') 
        self.title = name if title == None else title
        self.datapath = self.options.indir # usually /eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/logs_obelix_setup/
        self.temperature = temperature # temperature for measurement, usually -20 C
        self.structures = list(filter(lambda s: s != "GCD" or not noGCD, self.options.structures)) # usuallly ['MOShalf','MOS2000', 'GCD']
        self.configdir = self.options.configdir
        self.plotdir = f"{self.options.outdir}/{self.typeNameGoodString}_{self.name}"
        self.GCD_cuts = []
        self.GCD_maxV = 85.0 # might be passed from outside
        if self.configdir:
            self.readGCDcuts()

        self.subfolderPerDose = {}
        self.doses = {s: self.getAllDoses(s) for s in self.structures} # some doses might have multiple measurements
        #self.graphs = dict.fromkeys(self.structures, None) # will have a graph per structure
        self.graphsVsDose = {s : None for s in self.structures} # will have a graph per structure, with Vfb or GCD current vs dose
        self.graphsVsDose_alt = {s : None for s in self.structures} # will have a graph per structure, with oxide charge density or surface velocity vs dose
        self.graphsSingleDose = {s : {d: None for d in self.doses[s]} for s in self.structures}
        self.MOScurrentSingleDose = {}
        if "EPI" in self.typeName:
            self.MOScurrentSingleDose = {s : {d: None for d in self.doses[s]} for s in self.structures if "MOS" in s}
        self.Cox = {s : None for s in self.structures if "MOS" in s} # to be filled later once graphs exist
        self.plotStyle = {"lc" : ROOT.kBlack, #line color
                          "lw" : 2,           # line width                          
                          "ls" : 1,           # line style
                          "ms" : 20,          # marker style 
                      }
        self.readData()

    def run(self):
        createPlotDirAndCopyPhp(self.plotdir)
        self.plotData()
        self.makeAltGraph()
        #if "EPI" in self.typeName:
        #    self.plotMOScurrent() # when I run this function the final summary plot with all graphs becomes white, keep commented for now until I figure out what happens

    def getName(self):
        return self.name

    def getTypeName(self):
        #return getSampleTypeFromName(self.name)
        tmp = sast.getSampleAttribute(self.name, "leg")
        if hasattr(self.options, "samples") and self.options.samples != ["CandDose"] and "C EPI 14 kGy/h" in tmp:
            return "C EPI"
        else:
            return sast.getSampleAttribute(self.name, "leg")

    def getPlotStyle(self, key):
        return self.plotStyle[key]

    def setPlotStyle(self, key, val):
        self.plotStyle[key] = val

    def getGraphVsDose(self, structure, alt=False):
        return self.graphsVsDose_alt[structure] if alt else self.graphsVsDose[structure]

    def getGraphSingleDose(self, structure, dose):
        return self.graphsSingleDose[structure][dose]

    def getPath(self, mainpath, sample, structure, dose, temperatureFloat):
        trueDose = self.subfolderPerDose[dose][0]
        subfolderName = self.subfolderPerDose[dose][1]
        temperature = str(temperatureFloat).replace('-','m')
        prefix = "cv" if "MOS" in structure else "iv" if "GCD" in structure else "unknown" 
        path = f"{mainpath}/{sample}_{temperature}C_{trueDose}kGy"
        #path = glob(f"{path}/*/")[-1]
        path += f"/{subfolderName}/"
        path += f"{prefix}_{sample}_{temperature}C_{trueDose}kGy_{structure}.dat"
        #print(f"{sample} - {dose} kGy: {path}")
        if not os.path.exists(path):
            print(f"Error: I did not find this file {path}")
            quit()
        return path

    def getAllDoses(self, structure):
        temperature = str(self.temperature).replace('-','m')
        regexp = f"{self.name}_{temperature}C_[0-9]+kGy"
        #print(regexp)
        regMatch = re.compile(regexp)
        folders = [f for f in os.listdir(self.datapath) if regMatch.match(f)]
        #print(folders)
        doses = []
        for f in folders:
            dose = str(f.split("_")[-1]).rstrip("kGy")
            subfolders = [subf for subf in os.listdir(self.datapath + "/" + f)]
            subfolders = sorted(subfolders)
            if len(subfolders) > 1:
                print(f"{self.name} {structure} has {len(subfolders)} folders for dose {dose}: {subfolders}")
                if self.options.addRatio:
                    print(f"I will only use the first one because I have to make ratios with other graphs")
                    subfolders = [subfolders[0]] # need only a single folder, let's use the first measurement
            for i in range(len(subfolders)):
                #print(f"{structure}: {f}")
                if self.options.maxDose > 0 and int(dose) > self.options.maxDose: 
                    continue
                doseForGraph = int(dose) + 0.01 * i
                self.subfolderPerDose[doseForGraph] = [dose, subfolders[i]] # fill with true dose and subfolder name 
                doses.append(doseForGraph) # add just 1 kGy not to overlap points
        doses = sorted(doses)
        # print(f"{structure}: {doses}")
        return doses

    def getAnnealingPaths(self, structure):
        #prefix = "cv" if "MOS" in structure else "iv" if "GCD" in structure else "unknown" 
        regexp = f"{self.name}_.+C_[0-9]+kGy_annealingStep[0-9]+"
        regMatch = re.compile(regexp)
        folders = [f for f in os.listdir(self.datapath) if regMatch.match(f)]
        annealingSteps = {} # step : timestamp(date, time)
        for f in folders:
            #tmpdose = int(str(f.split("kGy_")[0]).split("_")[-1])
            annStep = str(f.split("_")[-1]).lstrip("annealingStep")
            subfolders = str(os.listdir(self.datapath + "/" + f)[0])
            timestamp = subfolders.split("_")[1:]
            annealingSteps[annStep] = timestamp
        # print(f"{structure}: {doses}")
        return annealingSteps


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

        xH = min(self.GCD_maxV, xH) # better not to go above about self.GCD_maxV V

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
                        meas.append(1.e12*float(data[-3])) # capacitance in pF (usually order of 1e-10 F in the file)
                        err.append(0)
                        current.append(1.e06*float(data[-1])) # current in microA
                    elif 'GCD' in structure:
                        meas.append(-1.e09*float(data[2])) # in nA
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
                # do not fill graphs vs dose at dose = 0, this value is not interesting
                if dose < 1.0:
                    continue 
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
                        tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, deltaVfb)
                        #tgeVsDose.SetPointError(tgeVsDose.GetN()-1, 0, 0) # no error on Y for now
                else:
                    [current, err] = self.getGCDcurrent(dose)
                    # do not use the point for dose = 0 for GCD
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

    def makeAltGraph(self):
        for structure in self.structures: 
            self.graphsVsDose_alt[structure] = self.graphsVsDose[structure].Clone(f"{self.graphsVsDose[structure].GetName()}_alt")
            allXvalues = list(self.graphsVsDose[structure].GetX())
            allYvalues = list(self.graphsVsDose[structure].GetY())
            tmp = self.graphsVsDose_alt[structure]
            if "MOS" in structure:
                for i,x in enumerate(allXvalues):
                    Vfb = allYvalues[i]
                    Nox = getOxideChargeDensity(Vfb, structure, self.Cox[structure])
                    tmp.SetPoint(i, x, Nox)
                    tmp.SetPointError(i, self.graphsVsDose[structure].GetErrorX(i), self.graphsVsDose[structure].GetErrorY(i) * Nox / Vfb)
                tmp.GetYaxis().SetTitle("Oxide charge density [cm^{-2 }]")
            else:
                for i,x in enumerate(allXvalues):
                    current = allYvalues[i]
                    surfV = getSurfaceVelocity(current)
                    tmp.SetPoint(i, x, surfV)
                    tmp.SetPointError(i, self.graphsVsDose[structure].GetErrorX(i), self.graphsVsDose[structure].GetErrorY(i) * surfV / current)
                tmp.GetYaxis().SetTitle("Surface velocity [cm/s]")
        return

    def plotMOScurrent(self):

        c = ROOT.TCanvas("mosCurrentCanvas","")
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

        c.Clear()  # this apparently does very bad things, it deletes canvas made afterwards !
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

        plotnameTag = f"{self.plotdir}/fit_{self.name}_{structure}_{dose}kGy"
        slope = getDerivative(tge)
        slope.GetXaxis().SetTitle("Voltage [-V]")
        slope.GetYaxis().SetTitle("#Delta C / #Delta V [pF/V]")
        slope.GetXaxis().SetTitleSize(0.045)
        slope.GetYaxis().SetTitleSize(0.045)
        slope.GetYaxis().SetTitleOffset(0.99)
        slope.SetMarkerStyle(20)
        #slope.Draw("apl")
        #nameNoExt = f"{plotnameTag}_firstDerivative"
        #for ext in ["png"]:
        #    c.SaveAs(f"{nameNoExt}.{ext}")

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
        deltaVfb = 0.5 * abs(Vfb_alt - Vfb)
        Vfb_bestEstimate = 0.5 * (Vfb + Vfb_alt)

        tge.SetMarkerStyle(20)
        tge.SetMarkerSize(0.3)
        tge.GetXaxis().SetTitle("Voltage [V-]")
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
        plat_ext.SetParameter(0, self.Cox[structure])
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
        lat.DrawLatex(0.45, 0.2, f"dose = {dose} kGy")
        lat.DrawLatex(0.45, 0.15, "flat-band voltage = {:0.1f} +/- {:0.1f} V".format(Vfb_bestEstimate, deltaVfb))

        nameNoExt = f"{plotnameTag}"
        for ext in ["png", "pdf"]:
            c.SaveAs(f"{nameNoExt}.{ext}")
        return (Vfb_bestEstimate, deltaVfb)


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
                        default="/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray_bkup01Nov2021/logs_obelix_setup/", help="Input folder")
    parser.add_argument("-o", "--outdir", type=str, 
                        default="./allplots", help="Output folder for plots")
    parser.add_argument("-c", "--configdir", type=str, 
                        default=None, help="Folder with configuration files")
    parser.add_argument("--max-dose", dest="maxDose", type=int, 
                        default=-1, help="Maximum dose to use for all samples, negative means to use all")
    # samples
    parser.add_argument("-s", "--samples", nargs="+", default=['N4791-1_LR','N4790-1_UL','N4791-6_UL','N4790-13_LR','N4789-10_UL','N4790-1_LR','N4790-13_UL','N4791-6_LR','N4788-9_LR'], help="List of samples to be used")
    parser.add_argument("--xs", "--exclude-samples-GCD", dest="excludeSamplesGCD", nargs="*", default=['N4789-12_UL'], help="List of samples to be excluded for GCD")
    #parser.add_argument("--samples", nargs="+", default=['N4791-1_LR','N4790-1_UL','N4791-6_UL'], help="List of samples to be used")
    parser.add_argument("--structures", nargs="+", default=['MOShalf','MOS2000','GCD'], help="List of structures to use")
    parser.add_argument("--compare-old-sample", dest="compareOldSample", action="store_true", default=False, help="Add comparison to an older sample (graphs read from root files directly)")
    parser.add_argument("-r", "--add-ratio", dest="addRatio", type=str, 
                        default=None, help="Add ratio plots using the sample passed here at denominator")
    parser.add_argument("--remove-large-ratio", dest="removeLargeRatio", type=float, 
                        default=-1, help="If positive, when a ratio would be above this value the dot is removed from the ratio plot")
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    allSamples = sast.getSampleNames()
    samples = []
    # just a bunch of hardcoded sets of sensors to access them more easily, hopefully just a temporary solution until I find a more intelligent way
    if args.samples == ["FZ"]:
        samples = [x for x in allSamples if "FZ" in sast.getSampleAttribute(x, "leg")]
    elif args.samples == ["FZnor"]:
        samples = [x for x in allSamples if "FZ" in sast.getSampleAttribute(x, "leg") and "#" not in sast.getSampleAttribute(x, "leg")]
    elif args.samples == ["EPI"]:
        samples = [x for x in allSamples if "EPI" in sast.getSampleAttribute(x, "leg")]
    elif args.samples == ["EPInor"]:
        samples = [x for x in allSamples if "EPI" in sast.getSampleAttribute(x, "leg") and "#" not in sast.getSampleAttribute(x, "leg")]
    elif args.samples == ["FZandEPI"]:
        samples = list(filter(lambda x: any(y in sast.getSampleAttribute(x, "leg") for y in ["FZ", "EPI"]), allSamples))
    elif args.samples == ["FZandEPInor"]:
        samples = list(filter(lambda x: all(y not in sast.getSampleAttribute(x, "leg") for y in ["7 kGy/h", "39 kGy/h", "#"]) and any(y in sast.getSampleAttribute(x, "leg") for y in ["FZ", "EPI"]), allSamples))
    elif args.samples == ["FZandEPInorAndDose"]:
        samples = list(filter(lambda x: "#" not in sast.getSampleAttribute(x, "leg") and any(y in sast.getSampleAttribute(x, "leg") for y in ["FZ", "EPI"]), allSamples))
    elif args.samples == ["allC"]:
        samples = list(filter(lambda x: "C " in sast.getSampleAttribute(x, "leg"), allSamples))
    elif args.samples == ["allCnoDose"]:
        samples = list(filter(lambda x: "C " in sast.getSampleAttribute(x, "leg") and not any(y in sast.getSampleAttribute(x, "leg") for y in ["7 kGy", "39 kGy"]), allSamples))
    elif args.samples == ["CandDose"]:
        samples = list(filter(lambda x: "C EPI" == sast.getSampleAttribute(x, "leg") or "kGy" in sast.getSampleAttribute(x, "leg"), allSamples))
    else:
        sample = args.samples
    
    samplesGCD = list(filter(lambda a: a not in args.excludeSamplesGCD, samples)) # might need to make it per sample, usually it is only for the GCD

    print()
    print(f"Analysing these structures: {args.structures}")
    print(f"Analysing {len(samples)} samples: {samples}")
    if len(samplesGCD) < len(samples):
        print(f">>> Excluding these for GCD: {args.excludeSamplesGCD}")
    print()

    sampleDict = {}
    for sample in samples:        
        sampleDict[sample] = siliconSensorSample(sample, -20, options=args, noGCD=(sample not in samplesGCD))
        sampleDict[sample].run()
        sampleDict[sample].setPlotStyle("lc", sast.getSampleAttribute(sample, "lc"))
        sampleDict[sample].setPlotStyle("ms", sast.getSampleAttribute(sample, "ms"))
        sampleDict[sample].setPlotStyle("ls", sast.getSampleAttribute(sample, "ls"))
        sampleDict[sample].setPlotStyle("lw", sast.getSampleAttribute(sample, "lw"))
    
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
    canvas.SetBottomMargin(0.12)
    canvas.cd()

    # repeat twice to make graphs having Vfb (GCD current) or oxide charge density (surface velocity) on the y axis
    for ivar in range(2):
        altGraph = True if ivar == 0 else False
        for structure in args.structures:
            colors = []
            markers = []
            lineStyles = []
            lineWidths = []
            graphs = []
            legendEntries = []
            ratioGraphs = []
            for sample in samples:                
                if structure == "GCD" and sample not in samplesGCD: 
                    continue
                graphs.append(        sampleDict[sample].getGraphVsDose(structure, alt=altGraph) )
                legendEntries.append( sampleDict[sample].getTypeName() )
                colors.append(        sampleDict[sample].getPlotStyle("lc") )
                markers.append(       sampleDict[sample].getPlotStyle("ms") )
                lineStyles.append(    sampleDict[sample].getPlotStyle("ls") )
                lineWidths.append(    sampleDict[sample].getPlotStyle("lw") )
            #
            # get some graphs already available from external files
            # at some point I should make this inclusion a little less hardcoded
            if args.compareOldSample:
                sampleExt = "N0538_25_LR"
                tf = safeOpenFile(f"tmpRootFiles/dose_{sampleExt}_{structure}.root")
                objname = None
                if "MOS" in structure:
                    objname = "gNox" if altGraph else "gVfb"
                else:
                    objname = "gJ" if altGraph else "gI"
                graphExt = safeGetObject(tf, objname, detach=False)
                tf.Close()
                graphs.append(graphExt)
                legendEntries.append(getSampleTypeFromName(sampleExt))
                colors.append(ROOT.kCyan+1)
                markers.append(20)
                lineStyles.append(1)
                lineWidths.append(2)
            #
            #print(graphs)
            #print(legendEntries)
            if args.addRatio:
                den = sampleDict[args.addRatio].getGraphVsDose(structure, alt=altGraph).Clone("den")
                #print([den.GetPointX(i) for i in range(den.GetN())])
                #print([den.GetPointY(i) for i in range(den.GetN())])
                for gr in graphs:
                    ratioGraphs.append(ROOT.TGraphErrors())
                    rg = ratioGraphs[-1]
                    rg.SetName(f"{gr.GetName()}_ratio")
                    for i in range(den.GetN()):
                        if den.GetPointX(i) == 0.0: # for ratios can exclude the points at dose = 0, since the measured values (e.g. Vfb for MOS) only depends on the sensor production process
                            continue
                        xval = den.GetPointX(i)
                        if xval > gr.GetPointX(gr.GetN()-1):
                            continue
                        ratio = gr.Eval(xval, 0, "S") / den.GetPointY(i)
                        # do not add dots with weird ratios, might just be bad measurements (usually it happens for dose = 0, but we remove it anyway)
                        if args.removeLargeRatio < 0.0 or ratio < args.removeLargeRatio:
                            rg.SetPoint(i, xval, ratio)

            # ranges for standard graphs (minX should be the same as for standard graphs)
            maxX = 0
            maxY = 0
            minY = 0
            for igr,gr in enumerate(graphs):
                maxX = max(maxX, gr.GetPointX(gr.GetN()-1))
                allY = list(gr.GetY())
                maxY = max(allY) if igr == 0 else max(maxY, max(allY))
                minY = min(allY) if igr == 0 else min(minY, min(allY))
            maxY *= 1.1
            useLogX = True
            maxX *= 3.0
        
            legWidth = 0.5 if args.samples == ["CandDose"] else 0.3
            canvas.SetTitle(structure)
            minX = 0.8
            if "GCD" in structure:
                #yTitle = "GCD current [nA]" # read from graphs directly
                minyLeg = 0.3
                maxyLeg = minyLeg + 0.06 * len(graphs)
                minxleg = 0.2
                maxxleg = minxleg + legWidth
                legCoords = f"{minxleg},{minyLeg},{maxxleg},{maxyLeg}"
                useLogY = False
                #legCoords = f"0.65,{minyLeg},0.95,{maxyLeg}"
            else:
                #yTitle = "Flat-band voltage [-V]" # read from graphs directly
                minyLeg = 0.12
                maxyLeg = minyLeg + 0.06 * len(graphs)
                maxxleg = 0.9
                minxleg = maxxleg - legWidth
                legCoords = f"{minxleg},{minyLeg},{maxxleg},{maxyLeg}"
                useLogY = False
                minY = 0

            xTitle = graphs[0].GetXaxis().GetTitle()
            yTitle = graphs[0].GetYaxis().GetTitle()
            #print(f"minX = {minX};   maxX = {maxX};   minY = {minY};   maxY = {maxY}")
            postfix = "_alt" if altGraph else ""
            ROOT.TGaxis.SetExponentOffset(-0.06, 0.00, "y") # X and Y offset for Y axis
            drawGraphs(graphs, f"{xTitle}::{minX},{maxX}", f"{yTitle}::{minY},{maxY}", f"summaryVsDose_compareSamples_{structure}{postfix}", outdirSummary, 
                       legendEntries=legendEntries, legendCoords=legCoords,
                       vecColors=colors, vecMarkerStyle=markers, vecLineStyle=lineStyles, vecLineWidth=lineWidths,
                       passCanvas=canvas, moreText=f"Structure: {structure}", useLogX=useLogX, useLogY=useLogY)

            if len(ratioGraphs):
                # ranges for ratios if any
                maxXratio = 0
                maxYratio = 0
                minYratio = 0
                for igr,gr in enumerate(ratioGraphs):
                    allY = list(gr.GetY())
                    maxXratio = max(maxXratio, gr.GetPointX(gr.GetN()-1))
                    maxYratio = max(allY) if igr == 0 else max(maxYratio, max(allY))
                    minYratio = min(allY) if igr == 0 else min(minYratio, min(allY))
                maxXratio *= 3.0
                #maxYratio *= 2.0 if "GCD" in structure else 1.8
                #minYratio *= 0.9
                diffY = maxYratio - minYratio
                minYratio -= 0.5 * diffY
                minYratio = max(0.0, minYratio)
                maxYratio += 1.5 * diffY
                maxyLeg = 0.89
                minyLeg = maxyLeg - 0.06 * len(graphs)
                maxxleg = 0.9
                minxleg = maxxleg - legWidth
                legCoords = f"{minxleg},{minyLeg},{maxxleg},{maxyLeg}"
                yRatioTitle = str(yTitle.split("[")[0]) + f"ratio over {sampleDict[args.addRatio].getTypeName()}"
                if args.samples == ["CandDose"] and "reference" in sampleDict[args.addRatio].getTypeName():
                    yRatioTitle = str(yTitle.split("[")[0]) + f"ratio over reference"
                drawGraphs(ratioGraphs, f"{xTitle}::{minX},{maxXratio}", f"{yRatioTitle}::{minYratio},{maxYratio}", f"summaryVsDose_compareSamples_{structure}{postfix}_ratio", outdirSummary, 
                           legendEntries=legendEntries, legendCoords=legCoords,
                           vecColors=colors, vecMarkerStyle=markers, vecLineStyle=lineStyles, vecLineWidth=lineWidths,
                           passCanvas=canvas, moreText=f"Structure: {structure}", useLogX=useLogX, useLogY=useLogY)


            # repeat with log scale on Y axis too
            if "GCD" in structure:
                useLogY = True
                minyLeg = 0.12
                maxyLeg = minyLeg + 0.06 * len(graphs)
                maxxleg = 0.95
                minxleg = maxxleg - legWidth
                legCoords = f"{minxleg},{minyLeg},{maxxleg},{maxyLeg}"
                maxY *= 1.5 # 3.0
                minY = 0.8 * minY
                drawGraphs(graphs, f"{xTitle}::{minX},{maxX}", f"{yTitle}::{minY},{maxY}", f"summaryVsDose_compareSamples_{structure}{postfix}_logY", outdirSummary, 
                           legendEntries=legendEntries, legendCoords=legCoords,
                           vecColors=colors, vecMarkerStyle=markers, vecLineStyle=lineStyles, vecLineWidth=lineWidths, 
                           passCanvas=canvas, moreText=f"Structure: {structure}", useLogX=useLogX, useLogY=useLogY)

