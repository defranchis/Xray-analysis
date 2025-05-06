#!/usr/bin/env python

import re, sys, os, os.path, subprocess, json, ROOT, copy, math
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

from CMS_lumi import *
from array import array
import shutil

#########################################################################

def printLine(marker='-', repeat=30):
    print(marker*repeat)

#########################################################################

def checkHistInFile(h, hname, fname, message=""):
    if not h:
        print("Error {msg}: I couldn't find histogram {h} in file {f}".format(msg=message,h=hname,f=fname))
        quit()

def safeGetObject(fileObject, objectName, quitOnFail=True, silent=False, detach=True):
    obj = fileObject.Get(objectName)
    if obj == None:
        if not silent:
            print(f"Error getting {objectName} from file {fileObject.GetName()}")
        if quitOnFail:
            quit()
    else:
        if detach and hasattr(obj, "SetDirectory"):
            obj.SetDirectory(0)
    return obj
        
def safeOpenFile(fileName, quitOnFail=True, silent=False, mode="READ"):
    fileObject = ROOT.TFile.Open(fileName, mode)
    if not fileObject or fileObject.IsZombie():
        if not silent:
            print(f"Error when opening file {fileName}")
        if quitOnFail:
            quit()
        else:
            return None
    elif not fileObject.IsOpen():
        if not silent:
            print(f"File {fileName} was not opened")
        if quitOnFail:
            quit()
        else:
            return None
    else:
        return fileObject
            
def checkNullObj(obj, objName="object", quitOnFail=True):

    if obj == None:
        print(f"Error with {objName}: it was None.")
        if quitOnFail:
            quit()
        else:
            return True
    else:
        return False
            
#########################################################################

def addStringToEnd(name, matchToAdd, notAddIfEndswithMatch=False):
    if notAddIfEndswithMatch and name.endswith(matchToAdd):
        return name
    elif not name.endswith(matchToAdd):
        return name + matchToAdd

#########################################################################

def graph_derivative(g):
    der = ROOT.TGraph()
    for i in range(0,g.GetN()-1):
        der.SetPoint(i,(g.GetPointX(i+1)+g.GetPointX(i))/2.,(g.GetPointY(i+1)-g.GetPointY(i))/(g.GetPointX(i+1)-g.GetPointX(i)))
    return der

def getDerivative(tge):
    res = ROOT.TGraphErrors()
    res.SetName(f"{tge.GetName()}_derivative")
    for i in range(tge.GetN()-1): # skip last point N because below we need to evaluate point N+1
        res.SetPoint(i, 0.5*(tge.GetPointX(i+1)+tge.GetPointX(i)), (tge.GetPointY(i+1)-tge.GetPointY(i)) / (tge.GetPointX(i+1)-tge.GetPointX(i)))
    return res

def cutGraph(old):
    new = ROOT.TGraphErrors()
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
    return 0.5 * (X[i] + X[i-1])

#########################################################################

def getZaxisReasonableExtremesTH2(h,nSigma=3,minZtoUse=None,maxZtoUse=None):

    htmp = ROOT.TH1D("htmp","",1000,h.GetMinimum(),h.GetMaximum())
    nbins = h.GetNbinsX() * h.GetNbinsY()    
    for ibin in range (1,nbins+1):
        val = h.GetBinContent(ibin)
        canFill = True
        if minZtoUse != None:
            if val < minZtoUse: canFill = False
        if maxZtoUse != None:
            if val > maxZtoUse: canFill = False
        if canFill: htmp.Fill(val)

    mean = htmp.GetMean()
    stddev = htmp.GetStdDev()
    retmin = max(h.GetMinimum(),mean - nSigma*stddev)
    retmax = min(h.GetMaximum(),mean + nSigma*stddev)
    return retmin,retmax

#########################################################################


#########################################################################

def getMinMaxHisto(h, excludeEmpty=True, sumError=True, 
                   excludeUnderflow=True, excludeOverflow=True,
                   excludeMin=None, excludeMax=None):
    
    # Warning, fix this function, GetBinContent with TH2 is not that simple, there are the underflow and overflow in each row and column
    # must check whether bin is underflow or overflow
    # therefore, the global bin is obtained as the number of bins +2, multiplied for each axis

    # excludeEmpty = True exclude bins with content 0.0. Useful when a histogram is filled with values in, for example, [1,2] but hassome empty bins
    # excludeMin/Max are used to select a range in which to look for maximum and minimum, useful to reject outliers, crazy or empty bins and so on
    # for histograms with non-negative values, excludeEmpty=True is equal to excludeMin==0.0

    # sumError is used to add or subtract error when looking for min/max (to have full error band in range)
    # when using excludeMin/Max, the errors are still ignored when evaluating the range

    # the better combination of options depends on dimension: for a TH1 is useful to visualize the error band in the plot range, while for a TH2 
    # only the bin content is interesting in the plot (the error is not reported with TH2::Draw, unless plotting it in a 3D space

    # one might exploit excludeMin/Max to select a rage depending on the distribution on the histogram bin content
    # for example, one can pass excludeMin=h.GetMean()-2*h.GetStdDev() and excludeMax=h.GetMean()+2*h.GetStdDev() so to 
    # select a range of 2 sigma around the mean

    dim = h.GetDimension()
    nbins = 0
    if   dim == 1: nbins = h.GetNbinsX() + 2
    elif dim == 2: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2)
    elif dim == 3: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2) * (h.GetNbinsZ() + 2)
    else:
        logging.error("In getMaxHisto(): dim = %d is not supported. Exit" % dim)
        quit()

    maxval = -sys.float_info.max
    minval = sys.float_info.max
    firstValidBin = -1
    for ibin in range (1,nbins+1):
        if excludeUnderflow and h.IsBinUnderflow(ibin): continue
        if excludeOverflow and h.IsBinOverflow(ibin): continue
        tmpmax = h.GetBinContent(ibin)
        tmpmin = h.GetBinContent(ibin)
        if excludeEmpty and tmpmin == 0.0: continue
        if excludeMin != None and tmpmin <= excludeMin: continue
        if excludeMax != None and tmpmax >= excludeMax: continue
        if firstValidBin < 0: 
            logging.debug("ibin %d:   tmpmin,tmpmax = %.2f, %.2f" % (ibin,tmpmin,tmpmax))
            firstValidBin = ibin
        if sumError:
            tmpmin -= h.GetBinError(ibin)
            tmpmax += h.GetBinError(ibin)
        if firstValidBin > 0 and ibin == firstValidBin:
            #the first time we pick a non empty bin, we set min and max to the histogram content in that bin
            minval = tmpmin
            maxval = tmpmax
            logging.debug("#### ibin %d:   min,max = %.2f, %.2f" % (ibin,minval,maxval))
        else:
            minval = min(minval,tmpmin)
            maxval = max(maxval,tmpmax)
        logging.debug("ibin %d:   min,max = %.2f, %.2f" % (ibin,minval,maxval))
    
    return minval,maxval

#########################################################################

def getMinimumTH(h, excludeMin=None):
    # get minimum excluding some values. For example, if an histogram has an empty bin, one might want to get the minimum such that it is > 0
    # underflow are not considered
    
    dim = h.GetDimension()
    retmin = sys.float_info.max

    if dim == 1:
        for ix in range(1,h.GetNbinsX()+1):
            if retmin > h.GetBinContent(ix):
                if excludeMin != None:
                    if h.GetBinContent(ix) > excludeMin: retmin = h.GetBinContent(ix)
                else:
                    retmin = h.GetBinContent(ix)

    elif dim == 2:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                if retmin > h.GetBinContent(ix,iy):
                    if excludeMin != None:
                        if h.GetBinContent(ix,iy) > excludeMin: retmin = h.GetBinContent(ix,iy)
                    else:
                        retmin = h.GetBinContent(ix,iy)

    elif dim == 3:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                for iz in range(1,h.GetNbinsZ()+1):
                    if retmin > h.GetBinContent(ix,iy,iz):
                        if excludeMin != None:
                            if h.GetBinContent(ix,iy,iz) > excludeMin: retmin = h.GetBinContent(ix,iy,iz)
                        else:
                            retmin = h.GetBinContent(ix,iy,iz)
                            

    else:
        raise RuntimeError("Error in getMinimumTH(): unsupported histogram's dimension (%d)" % dim)

    return retmin

#########################################################################

def getMaximumTH(h, excludeMax=None):
    # get maximum excluding some values. For example, if an histogram has a crazy bin, one might want to get the maximum value that is lower than that
    # overflow are not considered
    
    dim = h.GetDimension()
    retmax = sys.float_info.min

    if dim == 1:
        for ix in range(1,h.GetNbinsX()+1):
            if retmax < h.GetBinContent(ix):
                if excludeMax != None:
                    if h.GetBinContent(ix) < excludeMax: retmax = h.GetBinContent(ix)
                else:
                    retmax = h.GetBinContent(ix)

    elif dim == 2:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                if retmax < h.GetBinContent(ix,iy):
                    if excludeMax != None:
                        if h.GetBinContent(ix,iy) < excludeMax: retmax = h.GetBinContent(ix,iy)                        
                    else:
                        retmax = h.GetBinContent(ix,iy)

    elif dim == 3:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                for iz in range(1,h.GetNbinsZ()+1):
                    if retmax < h.GetBinContent(ix,iy,iz):
                        if excludeMax != None:
                            if h.GetBinContent(ix,iy,iz) < excludeMax: retmax = h.GetBinContent(ix,iy,iz)    
                        else:
                            retmax = h.GetBinContent(ix,iy,iz)

    else:
        raise RuntimeError("Error in getMaximumTH(): unsupported histogram's dimension (%d)" % dim)

    return retmax

#########################################################################

def fillTH2fromTH2part(h2out, h2in,
                       xbinLow=1, ybinLow=1,
                       xbinHigh=None, ybinHigh=None,
                       xoffset=0, yoffset=0):

    if xbinHigh == None:
        xbinHigh = h2out.GetNbinsX()
    if ybinHigh == None:
        ybinHigh = h2out.GetNbinsY()

    for ix in range(xbinLow, 1 + xbinHigh):
        for iy in range(ybinLow, 1 + ybinHigh):
            content = h2in.GetBinContent(ix + xoffset, iy + yoffset)
            error   = h2in.GetBinError(  ix + xoffset, iy + yoffset)
            h2out.SetBinContent(ix, iy, content)
            h2out.SetBinError(  ix, iy, error)

#########################################################################

def fillTH3fromTH3part(h3out, h3in,
                       xbinLow=1, ybinLow=1, zbinLow=1,
                       xbinHigh=None, ybinHigh=None, zbinHigh=None,
                       xoffset=0, yoffset=0, zoffset=0):

    if xbinHigh == None:
        xbinHigh = h3out.GetNbinsX()
    if ybinHigh == None:
        ybinHigh = h3out.GetNbinsY()
    if zbinHigh == None:
        zbinHigh = h3out.GetNbinsZ()

    for ix in range(xbinLow, 1 + xbinHigh):
        for iy in range(ybinLow, 1 + ybinHigh):
            for iz in range(zbinLow, 1 + zbinHigh):
                content = h3in.GetBinContent(ix + xoffset, iy + yoffset, iz + zoffset)
                error   = h3in.GetBinError(  ix + xoffset, iy + yoffset, iz + zoffset)
                h3out.SetBinContent(ix, iy, iz, content)
                h3out.SetBinError(  ix, iy, iz, error)


#########################################################################

# can't this use TH3.Projection?
def fillTH2fromTH3zrange(h2, h3, zbinLow=1, zbinHigh=1):
    for ix in range(1, 1 + h2.GetNbinsX()):
        for iy in range(1, 1 + h2.GetNbinsY()):
            error = ROOT.Double(0)
            h2.SetBinContent(ix, iy, h3.IntegralAndError(ix, ix, iy, iy, zbinLow, zbinHigh, error))
            h2.SetBinError(ix,iy,error);

#########################################################################

def fillTH2fromTH3zbin(h2, h3, zbin=1):
    #fillTH2fromTH3zrange(h2, h3, zbinLow=zbin, zbinHigh=zbin)

    for ix in range(1, 1 + h2.GetNbinsX()):
        for iy in range(1, 1 + h2.GetNbinsY()):
            h2.SetBinContent(ix, iy, h3.GetBinContent(ix, iy, zbin))
            h2.SetBinError(ix, iy, h3.GetBinError(ix, iy, zbin));

#########################################################################

def getTH2fromTH3(hist3D, name, binStart, binEnd=None):
    if binEnd == None:
        binEnd = binStart
    hist3D.GetZaxis().SetRange(binStart,binEnd)
    #print(f"getTH2fromTH3(): {name}   projecting bins from {binStart} to {binEnd}")
    # Order yx matters to have consistent axes!
    hist2D = hist3D.Project3D("yxe") # yxe is to make TH2 with y axis versus x axis 
    hist2D.SetName(name)
    return hist2D

    
#########################################################################

def fillTH3binFromTH2(h3, h2, zbin, scaleFactor=None):
    for ix in range(1, 1 + h2.GetNbinsX()):
        for iy in range(1, 1 + h2.GetNbinsY()):
            val   = h2.GetBinContent(ix, iy)
            error = h2.GetBinError(ix, iy)
            if scaleFactor != None:
                val   *= scaleFactor
                error *= scaleFactor
            h3.SetBinContent(ix, iy, zbin, val)
            h3.SetBinError(ix, iy, zbin, error);
            
#########################################################################

def createPlotDirAndCopyPhp(outdir, copyPHP=True, phpPath="./templates/index.php"):
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        if copyPHP:
            shutil.copy(phpPath, outdir)


#########################################################################

def getAxisRangeFromUser(axisNameTmp="", 
                         separator="::", 
                         rangeSeparator=","
                         ):
  
    setXAxisRangeFromUser = False;
    fields = axisNameTmp.split(separator)
    axisName = fields[0]
    
    if len(fields) > 1:
        setXAxisRangeFromUser = True;
        xmin = float(fields[1].split(rangeSeparator)[0])
        xmax = float(fields[1].split(rangeSeparator)[1])
    else:
        xmin = 0
        xmax = 0
        
    return axisName,setXAxisRangeFromUser,xmin,xmax


#########################################################################

def adjustSettings_CMS_lumi():

    ## dummy function to be called before using any other fucntion calling CMS_lumi
    ## for some reason, the settings of the very first plot are screwed up.
    ## To fix this issue, it is enough to call it to a dummy plot
    dummy = ROOT.TH1D("dummy","",10,0,10)
    for i in range(1,1+dummy.GetNbinsX()):
        dummy.SetBinContent(i,i)
    dummy.GetXaxis().SetTitle("x axis")
    dummy.GetYaxis().SetTitle("y axis")
    cdummy = ROOT.TCanvas("cdummy","",600,600)
    dummy.Draw("HE")
    CMS_lumi(cdummy,"",True,False)
    setTDRStyle()        
    ## no need to save the canvas    


#########################################################################

def drawTH1(htmp,
            labelXtmp="xaxis",
            labelYtmp="Events",
            outdir= "./",
            prefix = "distribution",
            outhistname = "",
            canvasSize="700,625",
            passCanvas=None,
            moreTextLatex="",
            skipTdrStyle=False
            ):



    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)

    h = htmp.Clone("htmp")

    cw,ch = canvasSize.split(',')
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetBottomMargin(0.14)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()

    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle(labelX)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetYaxis().SetTitle(labelY)
    h.GetYaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.04)
    if (setXAxisRangeFromUser): h.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h.GetYaxis().SetRangeUser(ymin,ymax)
    # force drawing stat box
    h.SetStats(1)
    h.Draw("HIST")
    canvas.RedrawAxis("sameaxis")
    if not skipTdrStyle: 
        setTDRStyle()
    # force drawing stat box
    ROOT.gStyle.SetOptStat(111110)
    ROOT.gStyle.SetOptFit(1102)
    #    
    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

    for ext in ["png","pdf"]:
        canvas.SaveAs(outdir + "{pfx}{hname}.{ext}".format(pfx=prefix,hname=("_"+outhistname) if len(outhistname) else "",ext=ext))



#########################################################################




# function to draw 2D histograms, can also plot profile along X on top
def drawCorrelationPlot(h2D_tmp,
                        labelXtmp="xaxis", labelYtmp="yaxis", labelZtmp="zaxis",
                        canvasName="default", plotLabel="", outdir="./",
                        rebinFactorX=0,
                        rebinFactorY=0,
                        smoothPlot=False,
                        drawProfileX=False,
                        scaleToUnitArea=False,
                        draw_both0_noLog1_onlyLog2=1,
                        leftMargin=0.16,
                        rightMargin=0.20,
                        nContours=51,
                        palette=55,
                        invertePalette=False,
                        canvasSize="700,625",
                        passCanvas=None,
                        bottomMargin=0.1,
                        plotError=False,
                        plotRelativeError=False,
                        lumi=None,
                        drawOption = "colz",
                        skipLumi=False
                        ):


    ROOT.TH1.SetDefaultSumw2()
    adjustSettings_CMS_lumi()

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h2D_tmp.RebinX(rebinFactorX)
        else:                             h2D_tmp.RebinX(len(rebinFactorX)-1,"",array('d',rebinFactorX)) # case in which rebinFactorX is a list of bin edges

    if (rebinFactorY): 
        if isinstance(rebinFactorY, int): h2D_tmp.RebinY(rebinFactorY)
        else:                             h2D_tmp.RebinY(len(rebinFactorY)-1,"",array('d',rebinFactorY)) # case in which rebinFactorX is a list of bin edges

    if plotError or plotRelativeError:
        herr = h2D_tmp.Clone(h2D_tmp.GetName()+"_err")
        herr.Reset("ICESM")
        for i in range(1,herr.GetNbinsX()+1):
            for j in range(1,herr.GetNbinsY()+1):
                errval = h2D_tmp.GetBinError(i,j)
                if plotRelativeError:
                    if h2D_tmp.GetBinContent(i,j) != 0.0:
                        errval = errval/h2D_tmp.GetBinContent(i,j)
                    else:
                        errval = 1.0 if errval == 0 else 0.0
                herr.SetBinContent(i,j,errval)
        h2D = herr
    else:
        h2D = h2D_tmp

    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),        
                                         ##array ("d", [0.70, 1.00, 0.34]),        
                                         ##array ("d", [0.00, 1.00, 0.82]),        
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)

    if palette > 0:
        ROOT.gStyle.SetPalette(palette)  # 55:raibow palette ; 57: kBird (blue to yellow, default) ; 107 kVisibleSpectrum ; 77 kDarkRainBow 
    ROOT.gStyle.SetNumberContours(nContours) # default is 20 
    if invertePalette:
        ROOT.TColor.InvertPalette()

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    labelZ,setZAxisRangeFromUser,zmin,zmax = getAxisRangeFromUser(labelZtmp)
        
    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)    
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetBottomMargin(bottomMargin)
    canvas.cd()

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)
    # normalize to 1
    if (scaleToUnitArea): h2D.Scale(1./h2D.Integral())

    h2DGraph = 0

    h2DPlot = 0
    if (not smoothPlot): h2DPlot = h2D
    else:
        h2DGraph = ROOT.TGraph2D()
        h2DGraph.SetNpx(300)
        h2DGraph.SetNpy(300)
        nPoint = 0
        for iBinX in range (1,1+h2D.GetNbinsX()):
            for iBinY in range(1,1+h2D.GetNbinsY()):
                h2DGraph.SetPoint(nPoint,h2D.GetXaxis().GetBinCenter(iBinX),h2D.GetYaxis().GetBinCenter(iBinY),h2D.GetBinContent(iBinX,iBinY))
                nPoint += 1
            

        h2DPlot = h2DGraph.GetHistogram()

    if plotLabel == "ForceTitle":
        h2DPlot.SetTitle(h2D_tmp.GetTitle())
  
    h2DPlot.GetXaxis().SetTitle(labelX)
    h2DPlot.GetYaxis().SetTitle(labelY)
    h2DPlot.GetXaxis().SetTitleSize(0.05)
    h2DPlot.GetXaxis().SetLabelSize(0.04)
    h2DPlot.GetXaxis().SetTitleOffset(0.95) # 1.1 goes outside sometimes, maybe depends on root version or canvas width
    h2DPlot.GetYaxis().SetTitleSize(0.05)
    h2DPlot.GetYaxis().SetLabelSize(0.04)
    h2DPlot.GetYaxis().SetTitleOffset(1.1)
    h2DPlot.GetZaxis().SetTitleSize(0.05)
    h2DPlot.GetZaxis().SetLabelSize(0.04)
    h2DPlot.GetZaxis().SetTitleOffset(1.2)

    h2DPlot.GetZaxis().SetTitle(labelZ) 
    h2DPlot.Draw(drawOption)
    
    if (setXAxisRangeFromUser): h2DPlot.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h2DPlot.GetYaxis().SetRangeUser(ymin,ymax)
    if (setZAxisRangeFromUser): h2DPlot.GetZaxis().SetRangeUser(zmin,zmax)


    h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.4)


    h2DProfile = 0
    if drawProfileX:
        h2DProfile = h2D.ProfileX("%s_pfx" %h2D.GetName())
        h2DProfile.SetMarkerColor(ROOT.kBlack)
        h2DProfile.SetMarkerStyle(20)
        h2DProfile.SetMarkerSize(1)
        h2DProfile.Draw("EPsame")
        
    # not yet implemented
    setTDRStyle()
    if not skipLumi and not plotLabel == "ForceTitle": 
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)

    if plotLabel == "ForceTitle":
        ROOT.gStyle.SetOptTitle(1)        

    #h2DPlot.GetZaxis().SetMaxDigits(1)  #for N>99, should use scientific notation, I'd like to make it work only with negative exponential but haven't succeeded yet
    # canvas.Modified()
    # canvas.Update()

    leg = ROOT.TLegend(0.39,0.75,0.89,0.95)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    if plotLabel not in ["", "ForceTitle"]: leg.AddEntry(0,plotLabel,"")
    if drawProfileX: leg.AddEntry(0,"Correlation = %.2f" % h2DPlot.GetCorrelationFactor(),"")
    leg.Draw("same")

    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 1):
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        
    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 2):
        canvas.SetLogz()
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_logZ.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        canvas.SetLogz(0)


##########################################################


def drawSingleTH1(h1,
                  labelXtmp="xaxis", labelYtmp="yaxis",
                  canvasName="default", outdir="./",
                  rebinFactorX=0,
                  draw_both0_noLog1_onlyLog2=1,
                  leftMargin=0.15,
                  rightMargin=0.04,
                  labelRatioTmp="Rel.Unc.::0.5,1.5",
                  drawStatBox=False,
                  legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                  canvasSize="600,700",  # use X,Y to pass X and Y size     
                  lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                  drawLineLowerPanel="luminosity uncertainty::0.025", # if not empty, draw band at 1+ number after ::, and add legend with title
                  passCanvas=None,
                  lumi=None,
                  drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                  textForLines=[],                       
                  moreText="",
                  moreTextLatex=""
                  ):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
        # case in which rebinFactorX is a list of bin edges
        else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(0)

    #ymax = max(ymax, max(h1.GetBinContent(i)+h1.GetBinError(i) for i in range(1,h1.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        ymin *= 0.9
        ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        if ymin < 0: ymin = 0

    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("HIST")
    h1err = h1.Clone("h1err")
    h1err.SetFillColor(ROOT.kRed+2)
    h1err.SetFillStyle(3001)  # 3001 is better than 3002 for pdf, while 3002 is perfect for png
    #h1err.SetFillStyle(3002)
    #h1err.SetFillStyle(3005)
    h1err.Draw("E2same")
    #h1.Draw("HIST same")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    leg.AddEntry(h1,"Value","L")
    leg.AddEntry(h1err,"Uncertainty","F")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        h1.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines): bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 30)

    if len(drawVertLines):
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = h1.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. if "#eta" in textForLines[0] else 6.
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = (1.1)*ymax/2.  
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)


  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        ratio = h1.Clone("ratio")
        den_noerr = h1.Clone("den_noerr")
        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        ratio.Divide(den_noerr)
        ratio.SetFillColor(ROOT.kGray+1)
        #den_noerr.SetFillColor(ROOT.kGray)
        frame.Draw()
        ratio.SetMarkerSize(0)
        ratio.SetMarkerStyle(0) # important to remove dots at y = 1
        ratio.Draw("E2same")

        line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineColor(ROOT.kRed)
        line.SetLineWidth(1)
        line.Draw("Lsame")

        if drawLineLowerPanel:
            legEntry,yline = drawLineLowerPanel.split('::')
            line2 = ROOT.TF1("horiz_line_2",str(1+float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line3 = ROOT.TF1("horiz_line_3",str(1-float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line2.SetLineColor(ROOT.kBlue)
            line2.SetLineWidth(1)
            line2.Draw("Lsame")
            line3.SetLineColor(ROOT.kBlue)
            line3.SetLineWidth(1)
            line3.Draw("Lsame")
            x1leg2 = 0.2 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.5 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.25 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.35 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.AddEntry(line2,legEntry,"L")
            leg2.Draw("same")

        
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)


################################################################

def drawNTH1(hists=[],
             legEntries=[],
             labelXtmp="xaxis", labelYtmp="yaxis",
             canvasName="default",
             outdir="./",
             rebinFactorX=0,
             draw_both0_noLog1_onlyLog2=0,                  
             leftMargin=0.15,
             rightMargin=0.04,
             labelRatioTmp="Rel.Unc.::0.5,1.5",
             drawStatBox=False,
             legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
             canvasSize="600,700",  # use X,Y to pass X and Y size     
             lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
             drawLineLowerPanel="", # if not empty, draw band at 1+ number after ::, and add legend with title
             passCanvas=None,
             lumi=None,
             drawLumiLatex=False,
             skipLumi=False,
             drawVertLines="", # coordinates of x where to print line
             textForLines=[],                       
             moreText="",
             moreTextLatex="",
             onlyLineColor=False,
             drawErrorAll=False # default draws error only on first histogram
):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    if len(hists) != len(legEntries):
        logging.warning("In drawNTH1: #(hists) != #(legEntries). Abort")
        quit()

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
        # case in which rebinFactorX is a list of bin edges
        else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    adjustSettings_CMS_lumi()
    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)
    

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)
    else:
        canvas.SetBottomMargin(0.15)


    h1 = hists[0]
    hnums = [hists[i] for i in range(1,len(hists))]
    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    #h1.SetMarkerSize(0)

    colors = [ROOT.kRed+2, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kAzure+2, ROOT.kMagenta, ROOT.kBlack]
    for ic,h in enumerate(hnums):
        # h.SetLineColor(colors[ic])
        # h.SetFillColor(colors[ic])
        # if ic==0: h.SetFillStyle(3004)   
        # if ic==2: h.SetFillStyle(3002)   
        # h.SetFillColor(colors[ic])
        # h.SetMarkerSize(0)
        h.SetLineColor(colors[ic])
        h.SetMarkerSize(0)
        if not onlyLineColor:
            h.SetFillColor(colors[ic])
            if ic==0: 
                h.SetFillStyle(3004)   
            if ic==1: 
                h.SetFillColor(0) 
                h.SetLineWidth(2) 
            if ic==2: 
                h.SetFillStyle(3002)           
            if ic==3:
                h.SetFillColor(0)
                h1.SetMarkerColor(ROOT.kGray+3)
                h1.SetMarkerStyle(25)
                #h1.SetMarkerSize(2)
        else:
            h.SetLineWidth(2)
    #ymax = max(ymax, max(h1.GetBinContent(i)+h1.GetBinError(i) for i in range(1,h1.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        # ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        # ymin *= 0.9
        # ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        # if ymin < 0: ymin = 0
        #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (h1.GetName(),ymin,ymax)
        ymin = 9999.9
        ymax = -9999.9
        for h in hists:
            if h.GetBinContent(h.GetMaximumBin()) > ymax: ymax = h.GetBinContent(h.GetMaximumBin())
            if h.GetBinContent(h.GetMinimumBin()) < ymin: ymin = h.GetBinContent(h.GetMinimumBin())
        if ymin < 0: ymin = 0
        ymax *= 1.2
        
    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("PE")
    for h in hnums:
        h.Draw("HE SAME" if drawErrorAll else "HIST SAME")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    for il,le in enumerate(legEntries):
        leg.AddEntry(hists[il],le,"PE" if il == 0 else "L" if onlyLineColor else "FL")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        for htmp in hists:
            htmp.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines): bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 30)

    if len(drawVertLines):
        nLines = len(drawVertLines)
        #offsetXaxisHist = h1.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 10.
        for i in range(nLines):
            #vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
            vertline.DrawLine(drawVertLines[i],0,drawVertLines[i],ymax)
        if len(textForLines):
            for i in range(len(textForLines)): # we need nLines
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = (1.1)*ymax/2.
                if i == 0:
                    bintext.DrawLatex(h1.GetXaxis().GetBinLowEdge(0) + sliceLabelOffset, ytext, textForLines[i])
                else:                    
                    bintext.DrawLatex(drawVertLines[i-1] + sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)


  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if not skipLumi:
        if not drawLumiLatex:
            if lumi != None: CMS_lumi(canvas,lumi,True,False)
            else:            CMS_lumi(canvas,"",True,False)
        else:
            latCMS = ROOT.TLatex()
            latCMS.SetNDC();
            latCMS.SetTextFont(42)
            latCMS.SetTextSize(0.045)
            latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
            if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
            else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)')

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        if len(hists) == 2:
            ratio = h1.Clone("ratio")
            den = hnums[0].Clone("den")
            den_noerr = hnums[0].Clone("den_noerr")
            for iBin in range (1,den_noerr.GetNbinsX()+1):
                den_noerr.SetBinError(iBin,0.)
            den.Divide(den_noerr)
            ratio.Divide(den_noerr)
            den.SetFillColor(ROOT.kGray)
            den.SetFillStyle(1001)
            #den_noerr.SetFillColor(ROOT.kGray)
            frame.Draw()
            frame.SetMarkerSize(0)
            frame.SetMarkerStyle(0) # important to remove dots at y = 1
            den.Draw("E2same")
            ratio.Draw("EPSAME")
        else:
            ratio = h1.Clone("ratio")
            den_noerr = h1.Clone("den_noerr")
            for iBin in range (1,den_noerr.GetNbinsX()+1):
                den_noerr.SetBinError(iBin,0.)

            ratio.Divide(den_noerr)
            ratio.SetFillColor(ROOT.kGray)
            ratio.SetFillStyle(1001)
            #den_noerr.SetFillColor(ROOT.kGray)
            frame.Draw()
            ratio.SetMarkerSize(0)
            ratio.SetMarkerStyle(0) # important to remove dots at y = 1
            ratio.Draw("E2same")

            ratios = []
            for i,h in enumerate(hnums):
                ratios.append(h.Clone("ratio_"+str(i+1)))
                ratios[-1].Divide(den_noerr)
                #ratios[-1].SetLineColor(h.GetLineColor())
                #ratios[-1].SetMarkerSize(0)
                #ratios[-1].SetMarkerStyle(0)
                #ratios[-1].SetFillColor(0)
                if h.GetFillColor():
                    ratios[-1].Draw("E2 SAME")
                else:
                    ratios[-1].Draw("HIST SAME")

        line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw("Lsame")

        if drawLineLowerPanel:
            legEntry,yline = drawLineLowerPanel.split('::')
            line2 = ROOT.TF1("horiz_line_2",str(1+float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line3 = ROOT.TF1("horiz_line_3",str(1-float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line2.SetLineColor(ROOT.kBlue)
            line2.SetLineWidth(1)
            line2.Draw("Lsame")
            line3.SetLineColor(ROOT.kBlue)
            line3.SetLineWidth(1)
            line3.Draw("Lsame")
            x1leg2 = 0.2 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.5 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.25 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.35 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.AddEntry(line2,legEntry,"L")
            leg2.Draw("same")

        
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)


################################################################

def drawGraphCMS(grList, 
                 xAxisNameTmp = "xAxis", 
		 yAxisNameTmp = "yAxis", 
		 canvasName = "default",
		 outputDIR = "./",
		 leg_roc = None, # text for legend
                 legendCoords = "0.5,0.15,0.9,0.35;2", # number after ; sets the number of columns
                 lumi = None,
                 vecMCcolors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kOrange+1, ROOT.kCyan+2, ROOT.kGray+2],
		 etabinText = "",
                 canvasSize="800,800",
                 passCanvas=None,
             ):


    adjustSettings_CMS_lumi()
    xAxisName = ""
    xmin = 0
    xmax = 0
    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(xAxisNameTmp)
    #
    yAxisName = ""
    ymin = 0
    ymax = 0
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(yAxisNameTmp)

    nGraphs = len(grList)

    cw,ch = canvasSize.split(',')
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetFillColor(0)
    canvas.SetGrid()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.06)
    canvas.cd()

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,ly1,lx2,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)

    for ig in range(0,nGraphs):
        grList[ig].SetMarkerStyle(20);
        grList[ig].SetMarkerColor(vecMCcolors[ig]);
        grList[ig].SetLineColor(vecMCcolors[ig]);
        grList[ig].SetLineWidth(2);
        grList[ig].SetFillColor(vecMCcolors[ig]);
        if ig == 0: 
            grList[ig].Draw("ap")
        else: 
            grList[ig].Draw("p same");
        leg.AddEntry(grList[ig],leg_roc[ig],"LF");        

    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    grList[0].GetXaxis().SetTitleSize(0.05);
    grList[0].GetXaxis().SetLabelSize(0.04);
    grList[0].GetYaxis().SetTitleOffset(1.3);
    grList[0].GetYaxis().SetTitleSize(0.05);
    grList[0].GetYaxis().SetLabelSize(0.04);
    grList[0].GetXaxis().SetTitle(xAxisName.c_str());
    grList[0].GetYaxis().SetTitle(yAxisName.c_str());
    if setXAxisRangeFromUser:
        grList[0].GetXaxis().SetRangeUser(xmin,xmax);
    if setYAxisRangeFromUser:
        grList[0].GetYaxis().SetRangeUser(ymin,ymax);

    setTDRStyle() # check if it doesn't screw things up
    if lumi != None: 
        CMS_lumi(canvas,lumi,True,False)
    else:   
        CMS_lumi(canvas,"",True,False)

    etabin = ROOT.TLatex()
    etabin.SetNDC(); # not sure it is needed
    etabin.SetTextSize(0.05);
    etabin.SetTextFont(42);
    etabin.SetTextColor(ROOT.kBlack);
    etabin.DrawLatex(0.15,0.15,etabinText);

    canvas.RedrawAxis("sameaxis");

    for ext in [".png",".pdf"]:
        canvas.SaveAs(outputDIR+canvasName+ext)


######


def drawGraphs(grList, 
               xAxisNameTmp = "xAxis", 
               yAxisNameTmp = "yAxis", 
               canvasName = "default",
               outputdir = "./",
               legendEntries = None, # text for legend
               legendCoords = "0.5,0.15,0.9,0.35;2", # number after ; sets the number of columns
               #lumi = None,
               vecColors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kOrange+1, ROOT.kCyan+2, ROOT.kGray+2],
               vecMarkerStyle = [],
               vecMarkerSize = [],
               vecLineStyle = [],
               vecLineWidth = [],
               canvasSize="800,800",
               passCanvas=None,
               moreText="",
               useLogX=False,
               useLogY=False
           ):


    adjustSettings_CMS_lumi()
    xAxisName = ""
    xmin = 0
    xmax = 0
    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(xAxisNameTmp)
    #
    yAxisName = ""
    ymin = 0
    ymax = 0
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(yAxisNameTmp)
    
    #print(f"minX = {xmin};   maxX = {xmax};   minY = {ymin};   maxY = {ymax}")
    frame = ROOT.TH1D("frame", moreText, 1, xmin, xmax)
    if useLogY:
        frame.SetBinContent(1,0.001)
    frame.SetMarkerSize(0)
    frame.SetMarkerColor(0)
    frame.SetLineColor(0)

    nGraphs = len(grList)
    markerStyle = vecMarkerStyle if len(vecMarkerStyle) else [20 for i in range(nGraphs)]
    markerSize = vecMarkerSize if len(vecMarkerSize) else [1.5 for i in range(nGraphs)]
    lineStyle = vecLineStyle if len(vecLineStyle) else [1 for i in range(nGraphs)]
    lineWidth = vecLineWidth if len(vecLineWidth) else [2 for i in range(nGraphs)]

    cw,ch = canvasSize.split(',')
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    if passCanvas == None:
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.cd()
        canvas.SetFillColor(0)
        canvas.SetGrid()
        canvas.SetLeftMargin(0.14)
        canvas.SetRightMargin(0.06)
        canvas.SetBottomMargin(0.12)
        canvas.cd()

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,ly1,lx2,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)

    frame.Draw()
    frame.SetStats(0)

    dictLegGraph = {} # to sort legend so that entries are in alphabetic order
    for ig in range(nGraphs):
        grList[ig].SetMarkerStyle(markerStyle[ig])
        grList[ig].SetMarkerSize(markerSize[ig])
        grList[ig].SetMarkerColor(vecColors[ig])
        grList[ig].SetLineColor(vecColors[ig])
        grList[ig].SetLineWidth(lineWidth[ig])
        grList[ig].SetLineStyle(lineStyle[ig])
        grList[ig].SetFillColor(vecColors[ig])
        #grList[ig].Draw("apl" if ig == 0 else "pl same")
        grList[ig].Draw("pl same")
        #print("Drawing graph !!!!!!!!")
        dictLegGraph[legendEntries[ig]] = grList[ig]
    for k in sorted(dictLegGraph.keys()):
        leg.AddEntry(dictLegGraph[k], k, "PL")       

     
    # better to set this outside of this function
    #ROOT.TGaxis.SetExponentOffset(-0.05, 0.00, "y") # X and Y offset for Y axis
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitle(xAxisName)
    frame.GetYaxis().SetTitle(yAxisName)
    if setXAxisRangeFromUser:
        frame.GetXaxis().SetRangeUser(xmin,xmax)
    if setYAxisRangeFromUser:
        frame.GetYaxis().SetRangeUser(ymin,ymax)

    canvas.RedrawAxis("sameaxis")
    leg.Draw("same")

    text = ROOT.TLatex()
    text.SetNDC() # not sure it is needed
    text.SetTextSize(0.05)
    text.SetTextFont(42)
    text.SetTextColor(ROOT.kBlack)
    text.DrawLatex(0.14, 0.92, moreText)

    if useLogX:
        canvas.SetLogx()
    else:
        canvas.SetLogx(0)

    if useLogY:
        canvas.SetLogy()
    else:
        canvas.SetLogy(0)

    for ext in ["png", "pdf"]:
        canvas.SaveAs(f"{outputdir}/{canvasName}.{ext}")

    # restore original version after saving plot
    if useLogX:
        canvas.SetLogx(0)
    if useLogY:
        canvas.SetLogy(0)

