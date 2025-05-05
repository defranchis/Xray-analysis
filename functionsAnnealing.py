## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from utility import cutGraph, findX
import os

def getSampleTypeFromName(name):
    sampleNames = {'N4791-1_LR'  : 'A FZ',
                   'N4790-1_UL'  : 'B FZ',
                   'N4790-1_LR'  : 'B FZ #2',
                   'N4791-6_UL'  : 'C FZ', 
                   'N4791-6_LR'  : 'C FZ #2',
                   'N4790-13_UL' : 'D FZ',
                   'N4790-13_LR' : 'D FZ #2',
                   'N4789-10_UL' : 'C EPI',
                   'N4788-9_LR'  : 'D EPI',
                   'N4789-10_LR' : 'C EPI #2',
                   'N4789-12_UL' : 'C EPI high dose',
                   'N4789-12_LR' : 'C EPI low dose',
                   'N0538_25_LR' : 'C FZ old'
    }
    return sampleNames[name]

        

# for annealing
def formPathAnnealing(generalInputDir, sample, structure, afterAnnElapsedTime,
                      annealingTime=60, # in min
                      annealingTemperature=60, # in C, as string
                      afterAnnTemperature=20, # in C, as string
                      extra=None):
    # structure = MOShalf or MOS2000 for now
    if 'MOS' not in structure:
        print(f"invalid structure: {structure}")
        return
    # e.g. 1011_LR_ann_60min_60C_MOShalf_20C_1h17min
    # or 1011_LR_ann_30min_60C_MOShalf_p20C_5h55min_1weekFreezer_1weekendRT
    name = f"{sample}_ann_{annealingTime}min_{annealingTemperature}C_{structure}_{afterAnnTemperature}C_{afterAnnElapsedTime}"
    if extra:
        name += f"_{extra}"
    indir = generalInputDir + name
    filenames = [f for f in os.listdir(indir) if os.path.isfile(os.path.join(indir, f)) and f.endswith(".cv")]
    if len(filenames) > 1:
        print(f"Warning: found more files in path {indir}:\n{filenames}\nPlease check!")
        quit()
    infile = indir + "/" + filenames[0]
    return infile

def getAllPathsWithMatch(generalInputDir, regexp):
    # e.g. 1011_LR_ann_60min_60C_MOShalf_20C_1h17min
    # or 1011_LR_ann_30min_60C_MOShalf_p20C_5h55min_1weekFreezer_1weekendRT
    regMatch = re.compile(regexp)
    folders = [f for f in os.listdir(generalInputDir) if regMatch.match(f)]
    if len(folders):
        return folders
    else:
        return None

def getAnnealingDetailsFromPath(string):
    sample,other = string.split("_ann_")
    tokens = other.split("_")
    annealingTime, annealingTemperature, structure, afterAnnTemperature  = tokens[:4]
    afterAnnElapsedTime  = tokens[4:]
    return sample,annealingTime,annealingTemperature,structure,afterAnnTemperature,afterAnnElapsedTime


# to read output file from annealing study
def readCV(filename):
    infile = open(filename, 'r')
    lines = infile.read().splitlines()
    fill = False
    V=[]
    C=[]
    for line in lines:
        if fill==True:
            if line == 'END': break
            info = line.split('\t')
            V.append(float(info[0]))
            C.append(float(info[1]))
        if line == 'BEGIN': fill = True
    return [V,C]

def annealingTimeFromStringToMinutes(stringTime):
    # examples 20min, 2h55min
    nh = 0
    nmin = 0
    if "h" in stringTime:
        nh,nmin = stringTime.split("h")
        nmin = nmin.replace("min","")
    else:
        nmin = stringTime.replace("min","")
    return int(nh) * 60 + int(nmin)

# def makePlotAnnealing(generalInputDir, sample, structure, afterAnnElapsedTime,
#                       annealingTime=60, # in min
#                       annealingTemperature=60, # in C, as string
#                       afterAnnTemperature=20, # in C, as string
#                       extra=None):

#     path = formPathAnnealing(generalInputDir, sample, structure, 
#                              afterAnnElapsedTime, annealingTime, 
#                              annealingTemperature, afterAnnTemperature, 
#                              extra)
def makePlotAnnealing(path, structure):

    V,C = readCV(path)
    tge = ROOT.TGraphErrors()
    for i in range(0, len(V)):
        Ci = C[i] * 1e12 # F to pF
        tge.SetPoint(i, V[i], Ci)
        tge.SetPointError(i, 0, 0)
    if 'MOS' in structure:
        tge.SetName('gC')
        tge.SetTitle('CV; V gate [V]; MOS C [pF]')
    elif 'GCD' in structure:
        tge.SetName('gI')
        tge.SetTitle('IV; V gate [V]; diode I [nA]')
    return tge

def fitVfbAnnealing(args, path, sample, structure, afterAnnTime):

    filenames = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.endswith(".cv")]
    if len(filenames) == 0:
        print(f"Warning: found no file in path {path}:\n{filenames}\nPlease check!")
        quit()
    elif len(filenames) > 1:
        print(f"Warning: found more files in path {path}:\n{filenames}\nPlease check!")
        quit()
    infile = path + "/" + filenames[0]
    tge = makePlotAnnealing(infile, structure)
    if tge == None: 
        return

    # here should process all folders, or get the Cox from outside
    # tge.Print()
    Cox = max(list(tge.GetY()))
    # print(f"Cox = {Cox}")
    
    cut = True
    if cut: 
        tge = cutGraph(tge)
    
    c = ROOT.TCanvas()

    minC = min(list(tge.GetY()))
    maxC = max(list(tge.GetY()))
    maxV = max(list(tge.GetX()))
    diff = maxC-minC

    #incomplete curves
    if structure =='MOS2000' and maxC<150: # was 200
        print(f"Doing MOS2000 but maximum capacity is {maxC} pF (<150, so incomplete curve), and min = {minC} pF")
        print("Will quit now")
        return

    low_ramp  = findX(minC + 0.5 * diff, tge)
    high_ramp = findX(minC + 0.9 * diff, tge)

    ramp = ROOT.TF1('ramp', 'pol1(0)', low_ramp, high_ramp)
    ramp.SetLineColor(ROOT.kRed+1)
    tge.Fit(ramp, 'q0+', '', low_ramp, high_ramp)

    low_ramp_alt  = findX(minC + 0.3 * diff, tge)
    high_ramp_alt = findX(minC + 0.95 * diff, tge)
    ramp_alt = ROOT.TF1('ramp_alt', 'pol1(0)', low_ramp_alt, high_ramp_alt) # may add + 10 to high_ramp_alt for better visualization, but still fit in range up to high_ramp_alt
    ramp_alt.SetLineColor(ROOT.kOrange+2)
    tge.Fit(ramp_alt, 'q0+', '', low_ramp_alt, high_ramp_alt)

    # Vfb = (ramp.GetParameter(0)-plat.GetParameter(0))/(ramp.GetParameter(1)-plat.GetParameter(1))*(-1.)
    Vfb = -1. * (ramp.GetParameter(0)-Cox) / ramp.GetParameter(1)
    Vfb_alt = -1. * (ramp_alt.GetParameter(0)-Cox) / ramp_alt.GetParameter(1)
    deltaVfb = abs(Vfb_alt - Vfb)

    tge.SetMarkerStyle(20)
    tge.SetMarkerSize(0.7)
    tge.GetXaxis().SetTitle("Voltage [V]")
    tge.GetYaxis().SetTitle("Capacitance [pF]")
    tge.Draw('ap')
    ramp_alt.Draw('l same')
    ramp.Draw('l same')
    
    ramp_ext = ROOT.TF1('ramp_ext', 'pol1(0)', high_ramp, Vfb+5)
    plat_ext = ROOT.TF1('plat_ext', 'pol0(0)', Vfb-10, maxV)
    ramp_ext.SetParameters(ramp.GetParameter(0), ramp.GetParameter(1))
    plat_ext.SetParameter(0, Cox)
    ramp_ext.SetLineColor(ROOT.kBlue)
    plat_ext.SetLineColor(ROOT.kBlue)

    ramp_ext.Draw('l same')
    plat_ext.Draw('l same')

    ## I tried to modify the linear function according to the parameters' uncertainty, 
    ## but it turns out the uncertainties are so small that the nominal and alternate curves basically overlaps
    ## so better to devise a different way to get an uncertainty on the measured Vfb
    # ramp_slopeUp = ROOT.TF1('ramp_slopeUp', 'pol1', low_ramp, high_ramp)
    # ramp_slopeUp.SetParameters(ramp.GetParameter(0) - ramp.GetParError(0),
    #                            ramp.GetParameter(1) + ramp.GetParError(1))
    # ramp_slopeUp.SetLineColor(ROOT.kOrange+2)
    # ramp_slopeUp.Draw('l same')

    l = ROOT.TLine(Vfb, minC, Vfb, maxC)
    l.SetLineColor(ROOT.kGreen+1)
    l.Draw('l same')

    lat = ROOT.TLatex()
    lat.SetNDC();
    lat.SetTextFont(42)
    lat.SetTextSize(0.04)
    lat.DrawLatex(0.54, 0.3, "flat-band voltage = {:0.1f} +/- {:0.1f}".format(Vfb, deltaVfb))
    
    nameNoExt = f"{args.outdir}/annealing_{sample}_{structure}_{afterAnnTime}"
    c.SaveAs(f"{nameNoExt}.png")
    # print(f"Vfb = {Vfb} +/- {deltaVfb}")
    return [Vfb, deltaVfb]


def processMOSannealing(args, sample, structure):

    if sample in args.MOS_exclude:
        return

    paths = getAllPathsWithMatch(args.indir, args.annealingPathRegexp)
    dict_afterAnnElapsedTime_Vfb = {}
 
    for p in paths:
        if sample not in p:
            continue
        if structure not in p:
            continue
        tokens = getAnnealingDetailsFromPath(p)
        annealingTime = tokens[1]
        annealingTemperature = tokens[2]
        afterAnnTemperature  = tokens[4]
        if len(tokens[5]) == 1: # for now remove other folders
            pass
            # print(p)
        else:
            continue
        afterAnnElapsedTime = tokens[5][0]
        ret = fitVfbAnnealing(args, args.indir+"/"+p, sample, structure, afterAnnElapsedTime)
        if ret == None:
            return
        else:
            Vfb = ret[0]
            deltaVfb = ret[1]
        
        afterAnnElapsedTimeMinutes = annealingTimeFromStringToMinutes(afterAnnElapsedTime)
        dict_afterAnnElapsedTime_Vfb[afterAnnElapsedTimeMinutes] = [Vfb, deltaVfb]

    allXvals = list(dict_afterAnnElapsedTime_Vfb.keys())
    if len(allXvals) == 0:
        return

    gVfb = ROOT.TGraphErrors()
    xmin = min(allXvals)
    xmax = max(allXvals)
    offsetX = 0.1 * (xmax - xmin)
    for itime,xtime in enumerate(sorted(allXvals)):
        # print((xtime, dict_afterAnnElapsedTime_Vfb[xtime]))
        gVfb.SetPoint(itime, xtime, dict_afterAnnElapsedTime_Vfb[xtime][0])
        gVfb.SetPointError(itime, 1, dict_afterAnnElapsedTime_Vfb[xtime][1])

    # gVfb.Print()
    vals = list(gVfb.GetY())        
    gmin = min(vals)
    gmax = max(vals)
    diff = gmax - gmin
    offsetY = (0.2 * diff) if structure == "MOS2000" else (0.1 * diff)
    if diff == 0.0:
        offsetY = 0.1 * gmax
    c = ROOT.TCanvas()
    gVfb.SetTitle('{} {}; time after annealing [min]; V flat-band [V]'.format(sample,structure))
    gVfb.SetMarkerStyle(7)
    # gVfb.GetYaxis().SetRangeUser(gmin - offsetY, gmax + offsetY) # this works but I do it before using the frame
    # gVfb.GetXaxis().SetRangeUser(xmin + offsetX, xmax + offsetX) # this doesn't seem to work as I would
    # draw TH1 to build frame and allow the user to extend the frame outside the graph range
    frame = ROOT.TH1F("frame","{} {}; time after annealing [min]; V flat-band [V]".format(sample,structure), 1, xmin - offsetX, xmax + offsetX)
    frame.GetXaxis().SetRangeUser(xmin - offsetX, xmax + offsetX)
    frame.GetYaxis().SetRangeUser(gmin - offsetY, gmax + offsetY)
    frame.SetStats(0)
    frame.Draw()
    gVfb.Draw('pl SAME')
    lat = ROOT.TLatex()
    lat.SetNDC();
    lat.SetTextFont(42)
    lat.SetTextSize(0.04)
    lat.DrawLatex(0.35, 0.3, f"{annealingTime} of annealing at {annealingTemperature}")
    lat.DrawLatex(0.35, 0.2, f"after annealing temperature = {afterAnnTemperature}")
    cName = f"{args.outdir}/VfbVersusTime_{sample}_{structure}"
    c.SaveAs(f"{cName}.png")

    return

# paths = getAllPathsWithMatch("/eos/user/h/hgsensor/HGCAL_test_results/Results_Xray/",".*_ann_30min_60C_MOShalf_20C.*")
# for p in paths:
#     tokens = getAnnealingDetailsFromPath(p)
#     if len(tokens[5]) == 1:
#         print(p)
# tokens = getAnnealingDetailsFromPath("1011_LR_ann_30min_60C_MOShalf_n20C_5h55min_1weekFreezer_1weekendRT_1weekRT_1nightRT")
# print(tokens)
# print(tokens[5][0])
