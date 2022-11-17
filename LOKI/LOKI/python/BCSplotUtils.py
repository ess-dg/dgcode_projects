from __future__ import print_function

#!/usr/bin/env python3
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
import os, matplotlib 
#matplotlib.use('agg')
##import PhysValTS.tns_style 
from PyAna import * 
#plt.enable_tex_fonts()

import numpy  as np

import SimpleHists as sh
from SimpleHists.plotutils import *
from SimpleHistsUtils.cmphists import cmphists as shu_cmphists
from matplotlib.colors import LogNorm
import PyAnaUtils.divide as divide
import sys, os
import math as m
from matplotlib.widgets import Slider
from matplotlib import gridspec
import re
import operator

#import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from cycler import cycler

from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import matplotlib.colors as colors
from matplotlib import colors as mcolors #added lated, might cause problems

from matplotlib import ticker
from matplotlib.pyplot import figure

from scipy.optimize import curve_fit
#Gaussian function for curve_fit
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def setPlottingParameters(multiplicationFactor = 1):
    #plt.rcParams["figure.figsize"] = (13,9.75)
    plt.rcParams['lines.linewidth'] = 4
    plt.rcParams['axes.linewidth'] = 1.5
    
    
    #plt.rcParams['font.size']= 10
    #plt.rcParams['axes.titlesize']= 2
    plt.enable_tex_fonts()

    SMALL_SIZE = 30 * multiplicationFactor #25 * multiplicationFactor
    LARGE_SIZE = 34 * multiplicationFactor #30 * multiplicationFactor
    plt.rc('font', size=LARGE_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=22)    # legend fontsize
    #plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    #plt.tight_layout()
    plt.rc('legend',**{'fontsize':20})


#from matplotlib import rc

#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

def _defcolmap(x):
    a, b = 0.15, 0.9
    return plt.cm.spectral(a + x * (b-a))

def _col(cmap, i, n):
    assert i < n
    return cmap(i / (n - 1.0)) if n > 1 else cmap(1.0)

cmap = _defcolmap


##########################################################################################################################################################
def shrinkArray(data, newDimX, newDimY, AxisX, AxisY):
    """This function rebins a numpy array (with summing values and get the mean of the axis values). If the new dimensions are not divisible with the old ones, 
    this function will APPEND enough ROWS AND/OR COLUMNS of ZEROS to make the dimensions divisible"""
    #print(newDimX, newDimY, 'new', AxisX.shape[0]/newDimX, AxisY.shape[0]/newDimY)
    #IF needed: fill the end of data and Axis arrays with zeros until the dimensions are divisible with newDimX and cols
    if (data.shape[0] % newDimX): 
        data = np.lib.pad(data, ((0,newDimX-(data.shape[0]%newDimX)), (0,0)), 'constant', constant_values = (0))
        AxisY = np.lib.pad(AxisY, ((0,newDimX-(AxisY.shape[0]%newDimX)), (0,0)), 'constant', constant_values = (AxisY[-1]))
    if (data.shape[1] % newDimY):
        data = np.lib.pad(data, ((0,0),(0,newDimY-(data.shape[1]%newDimY))), 'constant', constant_values=(0))
        AxisX = np.lib.pad(AxisX, ((0,0),(0,newDimY-(AxisX.shape[1]%newDimY))), 'constant', constant_values = (AxisX[-1]))

    data=data.reshape(newDimX, int(data.shape[0] / newDimX), newDimY, int(data.shape[1] / newDimY)).sum(axis = 1).sum(axis = 2) #sum the values along each axis
    AxisX=AxisX.reshape(newDimX, int(AxisX.shape[0]/newDimX), 1, 1).mean(axis = 1).sum(axis = 2) #get the mean of the axis values
    AxisY=AxisY.reshape(newDimY, int(AxisY.shape[0]/newDimY), 1, 1).mean(axis = 1).sum(axis = 2) 
    return [data, AxisX, AxisY]


##########################################################################################################################################################
def LegendSelector(shistFileNames, histName):
    """ Plot the same (1 dimensional) histogram from multiple shist files with the ability to hide or show any of them on the fly. Click on the line showing the color of a data in legend to show or hide that data plot"""

    plotNr = len(shistFileNames) #number of plots in the figure
    
    hc_dict = dict()
    histList = list()
    for i in range(plotNr):
        hc_dict[i] = sh.HistCollection(shistFileNames[i])
        histList.append(hc_dict[i].hist(histName))
   
    lArray=list()
    fig, ax = plt.subplots()

    
    for i in range(plotNr):
        histDataN, histErrorN, histAxisN, titleTextN, labelN = disassembleHist1d(histList[i])

        dXHalf = (histAxisN[1] - histAxisN[0]) / 2.0
        lIndex, = ax.plot(histAxisN, histDataN, lw=2, label = shistFileNames[i])
        
        #lIndex=ax.hist(frange(histAxisN[0]-dXHalf, histAxisN[-1]+dXHalf, 2*dXHalf), weights=histAxisN, label=shistFileNames[i])
        #lIndex=ax.hist(range(len(histAxisN)), weights=histAxisN, label=shistFileNames[i])
        lArray.append(lIndex)

        #plt.hist(range(2400, 2501), weights=a[0][0])
        
    ax.set_title(titleTextN)
    ax.set_xlabel(labelN)
    leg = ax.legend(fancybox = True, shadow = True)
        
    leg.get_frame().set_alpha(0.4)
    

    # we will set up a dict mapping legend line to orig line, and enable picking on the legend line
    lined = dict()
    for legline, origline in zip(leg.get_lines(), lArray):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

    def onpick(event):
        # on the pick event, find the orig line corresponding to the legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()
            
    fig.canvas.mpl_connect('pick_event', onpick)

   
    plt.show()
    
    return

##########################################################################################################################################################
def hist1DfromHist2D(hist2dim, index, isRow=False, normalise=False):
    """Make hist1D object from row or column of a hist2D object with the possibility of normlisation. Select row or column with index and isRow variables."""
    
    binContent2d = []
    binCenter2d = []
    
    if isRow:
        newHistTitle = hist2dim.getYLabel() + '=' + str(hist2dim.getBinCenterY(index))#+hist2dim.getTitle()
        h =  sh.Hist1D(newHistTitle, hist2dim.getNBinsX(), hist2dim.getXMin(), hist2dim.getXMax());
        h.setXLabel(hist2dim.getXLabel())
        for i in range(hist2dim.getNBinsX()):
            binContent2d.append(hist2dim.getBinContent(i, index))
            binCenter2d.append(hist2dim.getBinCenterX(i))
            
    else: #column
        newHistTitle = hist2dim.getXLabel() + '=' + str(hist2dim.getBinCenterX(index))#+hist2dim.getTitle()
        h = sh.Hist1D(newHistTitle, hist2dim.getNBinsY(), hist2dim.getYMin(), hist2dim.getYMax())
        h.setXLabel(hist2dim.getYLabel())
        for i in range(hist2dim.getNBinsY()):
            binContent2d.append(hist2dim.getBinContent(index, i))
            binCenter2d.append(hist2dim.getBinCenterY(i))

    if (normalise and sum(binContent2d)): # Normalise values if nonzero array (assuming no negative values)
        normFactor = (binCenter2d[1] - binCenter2d[0]) * sum(binContent2d)
        binContent2d = [x / normFactor for x in binContent2d]
            
            
    for j in range(len(binCenter2d)):
        h.fill(binCenter2d[j], binContent2d[j])
           
    h.setYLabel('Counts')
    h.setComment(hist2dim.getComment())
    #h.setErrorsByContent()

    return h

##########################################################################################################################################################
def hist2DFromData(histData, binCenterX, binCenterY, titleText, labelX, labelY, histComment=None):
    """Make a hist2D object from the values in a 2 dimensional array, the centers of the bins and labels along each axes and the title"""
    
    [nBinsX, nBinsY] = histData.shape
    
    Xmin = float(binCenterX[0])
    Xmax = float(binCenterX[-1])
    dx = float(binCenterX[1]) - float(binCenterX[0])
    Ymin = float(binCenterY[0])
    Ymax = float(binCenterY[-1])
    dy = float(binCenterY[1]) - float(binCenterY[0])
    
    h =  sh.Hist2D(nBinsX, Xmin - dx/2, Xmax + dx/2, nBinsY, Ymin - dy/2, Ymax + dy/2)
        
    for i in range(nBinsX):
        for j in range(nBinsY):
            h.fill(binCenterX[i], binCenterY[j], histData[i][j])

    h.setXLabel(labelX)
    h.setYLabel(labelY)
    h.setTitle(titleText)
    #h.setComment(histComment)

    return h

##########################################################################################################################################################
def overlayThemCmpHist(shistFileNames, histName, norm=1, logScale=1):
    """This function overlays arbitrary number of histograms from shist files"""
    hc_dict = dict()
    histograms_in = list()
    
    for i in range(len(shistFileNames)):
        hc_dict[i] = sh.HistCollection(shistFileNames[i])
        histograms_in.append((shistFileNames[i], hc_dict[i].hist(histName)))
        
    
    shu_cmphists(histograms_in,
                 rebin=None,
                 normalise=norm,
                 logy=logScale,)
   
    return

##########################################################################################################################################################
def overlayHistogramsInOneFileCmpHist(fileName, histNames, norm=1, logScale=1):
    """This function overlays arbitrary number of histograms from shist files"""
    hc = sh.HistCollection(fileName)
    histograms_in = list()
    
    for i in range(len(histNames)):
        histograms_in.append((histNames[i], hc.hist(histNames[i])))
        
    
    shu_cmphists(histograms_in,
                 rebin=None,
                 normalise=norm,
                 logy=logScale,)
   
    return


##########################################################################################################################################################
def overlayThem(shistFileNames, histName):
    """This function overlays arbitrary number of histograms from shist files"""
    
    hc_dict = dict()
    histList = list()
    for i in range(0, len(shistFileNames)):
        hc_dict[i] = sh.HistCollection(shistFileNames[i])
        histList.append(hc_dict[i].hist(histName))
                                                  
    legend=shistFileNames
    histList[0].overlay(histList[1 :]
                        , legend
                        , logy=True
                        , ymin=0
                        , title='Compare '+hc_dict[0].hist(histName).getTitle()
                        , xlabel=hc_dict[0].hist(histName).getXLabel()
                        , ylabel=hc_dict[0].hist(histName).getYLabel()
                        , errors=True
    )
    plt.show()
    return

##########################################################################################################################################################
def overlayThemHist1DList(histList, titleText=None, legend=None, logY=False):
    """This function overlays arbitrary number of hist1D histograms"""

    shistFileNames = []
    for i in range(len(histList)):
        shistFileNames.append(histList[i].getTitle())

    if not(titleText):  
        titleText = 'Compare '+histList[0].getTitle()
        

    #for i in range(len(shistFileNames)):
        #panelText, panelNumber = shistFileNames[i].split('=')
        ##numberOfNeutrons_MHz = str(int(histList[i].getBinSum(1,histList[i].getNBins())/1e6))
        ##shistFileNames[i] = "Panel "+str(int(float(panelNumber)+1))+" ("+numberOfNeutrons_MHz+" MHz)"
        #panelNumber=i
        #numberOfNeutrons_MHz = str( '%.2f' % float(histList[i].getBinSum(1,histList[i].getNBins())/1e6))
        #shistFileNames[i] = "Panel "+str(int(float(panelNumber)+1))+" ("+numberOfNeutrons_MHz+" MHz)"
        
    if not legend:
        legend=[shistFileNames[i] + '(Integral = '+ str(histList[i].getIntegral()) + ')' for i in range(len(histList))]

    #setPlottingParameters(multiplicationFactor = 1)
        
    softBracketXlabel = histList[0].getXLabel().replace('[','(')
    softBracketXlabel = softBracketXlabel.replace(']',')')
    softBracketXlabel = softBracketXlabel.replace('angstrom',r'\AA')
    softBracketXlabel = softBracketXlabel.replace('Wavelength',r'$\lambda$ ') 
    softBracketYlabel = histList[0].getYLabel().replace('[','(')
    softBracketYlabel = softBracketYlabel.replace(']',')')
    #softBracketYlabel = softBracketYlabel.replace('Counts','Neutrons/s')
    pdfName=titleText+'.pdf'
    titleText = '' #noTitle
    
    
    (fig1,ax1,leg1) = histList[0].overlay(histList[1 :]
                                      , legend
                                      , logy=logY 
                                      , ymin=0
                                      , title=titleText
                                      , xlabel=softBracketXlabel
                                      , ylabel=softBracketYlabel
                                      , errors=False
                                      , show=False
    )

    ax1.tick_params(direction='out', bottom=1, top=0, left=1, right=0, width=2, length=8 )

    
    #fig1, ax1 = plt.subplots()
    #print(pdfName)
    #plt.savefig(pdfName)
    
    #plt.savefig('myfig1.pdf', format='pdf', dpi=200)
    plt.show()
    
    return



def overlayConvDetEfficiency(filename, h_incidentName, h_conversionName, h_detectionName):

    hc = sh.HistCollection(filename)
    
    my_incident = hc.hist(h_incidentName)
    my_converted = hc.hist(h_conversionName)
    my_detected = hc.hist(h_detectionName) 
    
    h_conversion_efficiency = hist1dDivision(my_converted, my_incident, 'Conversion_efficiency')
    h_detection_efficiency = hist1dDivision(my_detected, my_incident, 'Detection_efficiency')

    overlayList = []
    overlayList.append(h_conversion_efficiency)
    overlayList.append(h_detection_efficiency)
    globalConversionEfficiency = float((my_converted.getIntegral() / my_incident.getIntegral()))
    globalDetectionEfficiency = float((my_detected.getIntegral() / my_incident.getIntegral()))
    conversionEfficiencyLegend = "Conversion efficiency (" + str( '%.2f' % (globalConversionEfficiency)) + ")"
    detectionEfficiencyLegend = "Detection efficiency (" + str( '%.2f' % (globalDetectionEfficiency)) + ")"
    
    overlayThemHist1DList(overlayList, titleText=None , legend=[conversionEfficiencyLegend, detectionEfficiencyLegend])

############################################################################################
def detectionToConversionRatio(filename, h_incidentName, h_conversionName, h_detectionName):

    hc = sh.HistCollection(filename)
    
    my_incident = hc.hist(h_incidentName)
    my_converted = hc.hist(h_conversionName)
    my_detected = hc.hist(h_detectionName) 

    h_detectionToConversionRatio = hist1dDivision(my_detected, my_converted, 'Detection To Conversion Ratio')

    nBins = h_detectionToConversionRatio.getNBins()
    
    histData = [h_detectionToConversionRatio.getBinContent(i) for i in range(nBins)]
    histAxis = [h_detectionToConversionRatio.getBinCenter(i) for i in range(nBins)]

    npData = np.array(histData)
    averageDTCRatio = np.average(npData[np.nonzero(npData)])

    npAxis=np.array(histAxis)

    histData = npData[np.nonzero(npData)]
    histAxis = npAxis[np.nonzero(npData)]
    

    
    setPlottingParameters(0.6)

    labelText = 'DCR values (' + str( '%.2f' % float((averageDTCRatio))) + ')'
    plt.plot(histAxis, histData, 'bo',  markersize=6, label=labelText)

    #labelTextAverage = 'Average (' + str( '%.2f' % float((averageDTCRatio))) + ')'
    #print labelTextAverage
    #plt.plot([0, 13], [averageDTCRatio, averageDTCRatio], 'r-', label=labelTextAverage)
    
    plt.legend( frameon=False)
    plt.ylabel('Detection to conversion ratio')
    plt.xlabel(r'$\lambda$(\AA)')
    plt.axis([0, 13, 0, 1])
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    plt.savefig('detectionToConversionRatio.pdf')
    
##########################################################################################################################################################
def hist2d4BasicOperations(hist1, operationSymbol, hist2):
    """Perform one of the 4 basic operations on two hist2d objects with matching bin centers. (addend1Hist, '+', addend2Hist) or (minuendHist, '-', subtrahendHist) or (multiplierHist, '*', multiplicandHist) or (dividendHist, '/', divisorHist) """
    
    hist1Data, hist1AxisX, hist1AxisY, hist1TitleText, hist1LabelX, hist1LabelY = disassembleHist2d(hist1)
    hist2Data, hist2AxisX, hist2AxisY, hist2TitleText, hist2LabelX, hist2LabelY = disassembleHist2d(hist2)
    
    if not (all(hist1AxisX == hist2AxisX)
            and all(hist1AxisY == hist2AxisY)):
        
        print("Error: hist1 and hist2 axes don't match!")
        return
    else:
        if(operationSymbol == '+'):
            operation = np.add
            operationText = 'Sum'
        elif(operationSymbol == '-'):
            operation = np.subtract
            operationText = 'Difference'
        elif(operationSymbol == '*'):
            operation = np.multiply
            operationText = 'Product'
        elif(operationSymbol == '/'):
            operation = np.divide
            operationText = 'Quotient'
            #indexes=np.where(hist2Data==0)
            
       
        resultData = operation(hist1Data, hist2Data)
        resultData = np.nan_to_num(resultData)#important after divison with zero
        resultTitleText = operationText+' of ('+hist1TitleText+' and '+hist2TitleText+')'
        
        resultHist = hist2DFromData(resultData, hist1AxisX, hist1AxisY, resultTitleText, hist1LabelX, hist1LabelY )
        
        #hist1.plot()
        #hist2.plot()
        #resultHist.plot()
        
        return resultHist

##########################################################################################################################################################
def disassembleHist2d(hist2d):
    """This function returns the Data, AxisX, AxisY, Title, LabelX and LabelY of a hist2D object in separate arrays/variables"""
    
    dimX = hist2d.getNBinsX()
    dimY = hist2d.getNBinsY()

    histData = np.empty([dimX, dimY], dtype=float)
    histAxisX = np.empty([dimX,1], dtype=float)
    histAxisY = np.empty([dimY,1], dtype=float)
    
    for j in range(dimY):
        for i in range(dimX):
            histData[i,j] = float(hist2d.getBinContent(i,j))
           
    for i in range(dimX):
        histAxisX[i] = float(hist2d.getBinCenterX(i))
    for i in range(dimY):
        histAxisY[i] = float(hist2d.getBinCenterY(i))


    labelX = hist2d.getXLabel()
    labelY = hist2d.getYLabel()
    titleText = hist2d.getTitle()
    
    return (histData, histAxisX, histAxisY, titleText, labelX, labelY) #comment?

##########################################################################################################################################################
def disassembleHist1d(hist1d):
    """This function returns the Data, Error, AxisX, Title and LabelX of a hist1D object in separate arrays/variables"""

    nBins = hist1d.getNBins()
    
    histData = [hist1d.getBinContent(i) for i in range(nBins)]
    histError = [hist1d.getBinError(i) for i in range(nBins)]
    histAxis = [hist1d.getBinCenter(i) for i in range(nBins)]
    titleText = hist1d.getTitle()
    label = hist1d.getXLabel()
    
    return (histData, histError, histAxis, titleText, label)

##########################################################################################################################################################
def OverlayAllRowOrCol(filename, histName, isRow=True, TitleText=None, newBinsX=None, newBinsY=None, plots=True, logY=False):
    """This function overlays all row or column of a hist2D object defined by hc file name and the histograms name. It is possible to set the title to override the default 'Compare *.getTitle()'. Shrinking the historam by summing rows and colums is also possible. """
    
    hc = sh.HistCollection(filename)
    hist2plot = hc.hist(histName)

    if(newBinsX and newBinsY):
        histData, histAxisX, histAxisY, hist2plotTitle, histLabelX, histLabelY = disassembleHist2d(hist2plot)
        histData, histAxisX, histAxisY = shrinkArray(histData, newBinsX, newBinsY, histAxisX, histAxisY)
        hist2plot = hist2DFromData(histData, histAxisX, histAxisY, hist2plot.getTitle(), histLabelX, histLabelY)

    binsX, binsY = hist2plot.getNBinsX(), hist2plot.getNBinsY()

    overlayList = []
    normOverlayList = []

    if isRow: #Overlay all Rows
        for i in range(binsY):
            overlayList.append( hist1DfromHist2D( hist2plot, i, isRow=True) ) # (newHist2plot,i,True)-->make hist1D from the i-th row of newHist2plot
            normOverlayList.append( hist1DfromHist2D( hist2plot, i, isRow=True, normalise='normalise') ) # normalise--> every hist1D object will be normalised to 1
        legendText=[str(i) for i in range(binsY)]
    else:
        for i in range(binsX):
            overlayList.append( hist1DfromHist2D( hist2plot, i, isRow=False) ) # (newHist2plot,i,False)-->make hist1D from the i-th column of newHist2plot
            normOverlayList.append( hist1DfromHist2D( hist2plot, i, isRow=False, normalise='normalise') ) # normalise--> every hist1D object will be normalised to 1
        legendText=[str(i) for i in range(binsX)]

    if(plots):
        #overlayThemHist1DList(normOverlayList,'Normalised ' + TitleText, legendText)
        overlayThemHist1DList(overlayList, TitleText , legendText, logY=logY)
        #overlayThemHist1DList(normOverlayList,TitleText)

        
    return overlayList
    #return normOverlayList

##########################################################################################################################################################
def OverlayDivisionAllRowOrCol(filename, histName1, histName2, isRow=True, TitleText=None, newBinsX=None, newBinsY=None, plots=True):
    """Temp function to plot Panel Efficiency"""

    hc = sh.HistCollection(filename)
    hist2plot1 = hc.hist(histName1)
    hist2plot2 = hc.hist(histName2)
    
    #if(newBinsX and newBinsY):
    #    histData, histAxisX, histAxisY, hist2plotTitle, histLabelX, histLabelY = disassembleHist2d(hist2plot)
    #    histData, histAxisX, histAxisY = shrinkArray(histData, newBinsX, newBinsY, histAxisX, histAxisY)
    #    hist2plot = hist2DFromData(histData, histAxisX, histAxisY, hist2plot.getTitle(), histLabelX, histLabelY)

    binsX, binsY = hist2plot1.getNBinsX(), hist2plot1.getNBinsY()

    overlayList1 = []
    overlayList2 = []


    if isRow: #Overlay all Rows
        for i in range(binsY):
            overlayList1.append( hist1DfromHist2D( hist2plot1, i, isRow=True) ) # (newHist2plot,i,True)-->make hist1D from the i-th row of newHist2plot
            overlayList2.append( hist1DfromHist2D( hist2plot2, i, isRow=True) ) # 
    else:
        for i in range(binsX):
            overlayList1.append( hist1DfromHist2D( hist2plot1, i, isRow=False) ) # (newHist2plot,i,False)-->make hist1D from the i-th column of newHist2plot
            overlayList2.append( hist1DfromHist2D( hist2plot2, i, isRow=False) ) 

    overlayList = []
    for i in range(5):
        overlayList.append(hist1dDivision(overlayList1[i], overlayList2[i], TitleText))

    if(plots):

        #overlayThemHist1DList(normOverlayList,'Normalised ' + TitleText)
        overlayThemHist1DList(overlayList, TitleText)
        ##overlayThemHist1DList(normOverlayList,TitleText)

    return overlayList
    #return normOverlayList


    
def safe_div(x,y):
    if y == 0:
        return 0
    elif x>y:
        return 1
    return x / y
    
##########################################################################################################################################################
def hist1dDivision(hist1d1, hist1d2, titleText):
    #"""This function returns the Data, Error, AxisX, Title and LabelX of a hist1D object in separate arrays/variables"""

    nBins = hist1d1.getNBins()

    
    
    histData = [safe_div(hist1d1.getBinContent(i),hist1d2.getBinContent(i)) for i in range(nBins)]
    histError = [hist1d1.getBinError(i) for i in range(nBins)]#NOT OK
    histAxis = [hist1d1.getBinCenter(i) for i in range(nBins)]
    #titleText = hist1d.getTitle()
    label = hist1d1.getXLabel()

    
    h =  sh.Hist1D(titleText, len(histAxis), hist1d1.getXMin(), hist1d1.getXMax());
    h.setXLabel(label)

                               
    for j in range(len(histAxis)):
        h.fill(histAxis[j], histData[j])
           
    h.setYLabel('Efficiency')
    #h.setComment(hist1d1.getComment())

    return h
    

##########################################################################################################################################################
def analyse2D(fileName, histName, newBinsX=None, newBinsY=None):
    """This function helps to analyse 2D plots showing a selected row and column of the histogram given by the fileName and the histName. By giving newnewBinsX and newBinsY values one can rebin(shrink) the histogram. Two sliders are provided to select the row or colunm to plot."""

    hc = sh.HistCollection(fileName)
    hist2plot = hc.hist(histName)
    
    histData, histAxisX, histAxisY, hist2plotTitle, histLabelX, histLabelY = disassembleHist2d(hist2plot)        

    if(newBinsX and newBinsY): #if user wants to rebin before ploting
        histData, histAxisX, histAxisY = shrinkArray(histData, newBinsX, newBinsY, histAxisX, histAxisY)
        hist2plot = hist2DFromData(histData, histAxisX, histAxisY, hist2plot.getTitle(), histLabelX, histLabelY)

    [dimX, dimY] = histData.shape
    histCenterY = dimY/2
    histCenterX = dimX/2
    
#----------------------------------------------
    
    fig = plt.figure() # set up figure
    gs = gridspec.GridSpec(2, 2)

    ax = fig.add_subplot(gs[0,:]) 
    hist2plot.plot(show=False, figure=fig, axes=ax) #plot the 2d data

    #----------------------------------------------
    bx = fig.add_subplot(223)
    bx.autoscale(True)
    
    # plot first data set (values at the center of the axes)
    frameX = histCenterY
    barPlot = bx.bar(histAxisX, histData[:,frameX], width = (histAxisX[2,0] - histAxisX[1,0]))
    bx.set_xlabel(histLabelX)
    bx.set_title(histLabelY + '=' + str(histAxisY[frameX]))
    
    # make the slider
    bxframe = plt.axes([0.15, 0.01, 0.3, 0.03])
    bframe = Slider(bxframe, histLabelY + ' step', 0, dimY, valinit = histCenterY, valfmt = '%d')
    
    # call back function
    def update(val):
        frameX = int(numpy.floor(bframe.val)) #index of the selected value

        for i in range(frameY):
            barPlot[i].set_height(histData[i,frameX]) #update the barPlot
            
        bx.set_title(histLabelY + '=' + str(histAxisY[frameX])) #update the title
        bx.relim()
        bx.autoscale_view()
        
       # bx.clear() #Not sure if it is faster or slower
       # bx.bar(histAxisX,histData[:,frameX],width=(histAxisX[2,0]-histAxisX[1,0]))
       # bx.set_title(histLabelY+'='+str(histAxisY[frameX]))
       # bx.set_xlabel(histLabelX)
        
    bframe.on_changed(update) # connect callback to slider
    
    #------------------------------------------------
    cx = fig.add_subplot(224)
    cx.autoscale(True)

    # plot first data set (values at the center of the axes)
    frameY = histCenterX 
    barPlot2 = cx.bar(histAxisY, histData[frameY,:], width = (histAxisY[2,0]-histAxisY[1,0]))
    cx.set_xlabel(histLabelY)
    cx.set_title(histLabelX + '=' + str(histAxisX[frameY]))
    
    # make the slider
    cxframe = plt.axes([0.6, 0.01, 0.3, 0.03])
    cframe = Slider(cxframe, histLabelX + ' step', 0, dimX, valinit = histCenterX, valfmt='%d')
    
    # call back function
    def update2(val):
        frameY = int(numpy.floor(cframe.val)) #index of the selected value
        
        for i in range(frameX):
            barPlot2[i].set_height(histData[frameY,i]) #update the barPlot
             
        cx.set_title(histLabelX + '=' + str(histAxisX[frameY])) #update the title
        cx.relim()
        cx.autoscale_view()
        
    cframe.on_changed(update2) # connect callback to slider  

    plt.show()
    return#---------------------------------------------


##########################################################################################################################################################
def visualizeRate(DummyHist,RateHist,titleText):

    DummyHistData, DummyHistAxisX, DummyHistAxisY, _, DummyLabelX, DummyLabelY = disassembleHist2d(DummyHist)
    convRateHistData = disassembleHist2d(RateHist)[0]

    [nBinsXDummy, nBinsYDummy] = DummyHistData.shape
    [nBinsX, nBinsY] = convRateHistData.shape
    
    for i in range(nBinsXDummy):
        for j in range(nBinsYDummy):
            if(DummyHistData[i,j] > 1):
                tubeNumber = int(DummyHistData[i,j] // 100)
                strawNumber = int(DummyHistData[i,j] % 100)
                
                DummyHistData[i,j] = convRateHistData[strawNumber / 10, tubeNumber]
            elif(DummyHistData[i,j] == 1): # exception because there is a straw where the tubCpNr and StrCpNr is zero and I made the Dummy with modified rata_analysis to set DummyHistData[i,j] one instead of zero for that particular straw
                DummyHistData[i,j]=convRateHistData[0,0]
                

    my_hist_visualizeRate = hist2DFromData(DummyHistData, DummyHistAxisX, DummyHistAxisY, titleText, DummyLabelX, DummyLabelY )
    

    return(my_hist_visualizeRate)

##########################################################################################################################################################
def newVisualizeRate(RateHist, titleText, logaritmicalPlot=1, numberOfTubes = 40, numberOfPanels = 5, xLabelText='y (mm)', yLabelText='z (mm)', colorbarLabelText='Neutrons/straw/ms', plotFull = 1, createPDF=0):
    #This function visualises the data from a rate histogram on the cross-section image of the bcs detector model with arbitrary number panels
    # And no more than 10 (currently, but upgrading to 100) tubes
    
    convRateHistData = disassembleHist2d(RateHist)[0]

    convRateHistData[convRateHistData>1e10]=1
    if titleText=='Efficiency_Geantino_Detection': #FOR GEANTINOS
        print('replace')
        convRateHistData[convRateHistData==0]=1  
    ###convRateHistData[convRateHistData != 1e10]=0.15 ###demonstrationFigure
    #convRateHistData[convRateHistData>0.5]=0

    print(convRateHistData)
    
    patches = []
    colors_files=[]

    numberOfStraws = 7
    straw_radius = 3.75; #[mm]

    #position of straw centers in a tube before making a phisical volume van a ratation applied. [mm]
    strawCenterBeforeRotation = np.array([[-2 * straw_radius, 0],
                                          [-straw_radius, -straw_radius * m.tan(m.pi/3)],
                                          [-straw_radius, straw_radius * m.tan(m.pi/3)],
                                          [0, 0],
                                          [straw_radius, -straw_radius * m.tan(m.pi/3)],
                                          [straw_radius, straw_radius * m.tan(m.pi/3)],
                                          [2 * straw_radius, 0]])

    #apply a rotation of alpha degrees on the center coordinates
    strawCenter = np.array([[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]])
    alpha = 50 * (m.pi/180) # 50 degree in radian 
    for i in range(7):
        strawCenter[i,0] = strawCenterBeforeRotation[i,0] * 1.0*m.cos(alpha) - strawCenterBeforeRotation[i,1] * m.sin(alpha) #
        strawCenter[i,1] = strawCenterBeforeRotation[i,0] * m.sin(alpha) + strawCenterBeforeRotation[i,1] * m.cos(alpha) #


    sample_detector_distance = 5012.7 # 5012.7 #[mm]
    tube_box = 25.4 # [mm] 
    panel_offset = 10.16 #[mm]
    
    sum_offset = 0 #shift of the adjacent panels
    maximumValue = 0
    averageValue = 0
    
    highestPanelSum = 0
    for panelNumber in range(numberOfPanels):
        panelSum=0
        for tubeNumber in range(numberOfTubes):

            if(sum_offset >= tube_box): # the offset cannot be larger than a tube_box 
                sum_offset -= tube_box  # we want the panels to cover each other

            yTubeCoord = (tubeNumber * tube_box - (numberOfTubes - 1) * 0.5 * tube_box + sum_offset) #coordinates of the tube center
            zTubeCoord = sample_detector_distance + panelNumber * tube_box
                
            for strawNumber in range(numberOfStraws):

                # there was a change in storing data to histograms when moving to 40 tubes instead of 7
                # the tube number [0,40[ and straw number goes from  [0,7[ so the unique ID is (tube number)*100 + straw number
                # instead of the old (tube number)*10 + straw number
                if(numberOfTubes<=10): 
                    tubeID = panelNumber * 10 + tubeNumber
                else:
                    tubeID = panelNumber * 100 + tubeNumber
                    
                #fill the matching rate values and the circle objects(representing straws) into the color_files and patches arrays
                colors_files.append(convRateHistData[strawNumber, tubeID])
                panelSum = panelSum+convRateHistData[strawNumber, tubeID];
                #colors_files.append(strawNumber) #debug
                 
                yCoord = yTubeCoord + strawCenter[strawNumber, 0]
                zCoord = zTubeCoord + strawCenter[strawNumber, 1]
                
                circle = Circle((yCoord,zCoord), straw_radius) #Circle((x1,y1), r)
                patches.append(circle)

                averageValue = averageValue + convRateHistData[strawNumber, tubeID]
                if (convRateHistData[strawNumber, tubeID]>maximumValue):
                    maximumValue = convRateHistData[strawNumber, tubeID]
                    #print ( str(maximumValue) + "  " + (str)(panelNumber+20) +'-'+ (str)(panelNumber+21)+" ms" )

        sum_offset += panel_offset
        
        if panelSum > highestPanelSum:
            highestPanelSum = panelSum
            #print(str(highestPanelSum)+" 'Panel' Number: "+str(panelNumber))

    averageValue = averageValue / (numberOfPanels * numberOfTubes * numberOfStraws) 
    #print("Global peak incident (inc/det) rate (from FirstPanel): " +str(highestPanelSum))

    if plotFull:
        setPlottingParameters(multiplicationFactor=1.8)
    else:
        setPlottingParameters()
        
    fig = plt.figure()
    ax = plt.gca()
    

    if(logaritmicalPlot):
        p_files = PatchCollection(patches, cmap=matplotlib.cm.jet, norm=colors.LogNorm(),  alpha=0.7)
    else:
        p_files = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.7) 

    p_files.set_array(np.array(colors_files))
    ax.add_collection(p_files)
    ax.set_title('')
    #ax.set_title(titleText)
    ax.axis('image')
    ax.set_ylabel(yLabelText.replace('(','[').replace(')',']'))
    ax.set_xlabel(xLabelText.replace('(','[').replace(')',']'))
    #ax.set_xlim([- (numberOfTubes/2*tube_box + 20), (numberOfTubes/2*tube_box + 35)])
    #ax.set_xlim([- (numberOfTubes/2*tube_box +5), (numberOfTubes/2*tube_box + 25)])


    if plotFull:
        ax.set_ylim([ 5000 - 5, 5000 +(tube_box * numberOfPanels) + 5]) #current best
        ax.set_xlim([-510, 530])
        clbar = plt.colorbar(p_files, orientation="horizontal", aspect=100, pad=0.20, fraction=0.1)
        fig.set_size_inches(45,10)#current best
    else:
        clbar = plt.colorbar(p_files, orientation="horizontal", aspect=50) ###OFF for demonstration figure
        ax.set_xlim([-125, 125])
        ##fig.set_size_inches(12,9)
        fig.set_size_inches(2*6.3889,2*4.7917)  ###OFF for demonstration figure
        ###fig.set_size_inches(2*6.3889,2*3.9)  ### for demonstration figure

    p_files.set_clim([0, 0.5])
    #tick_locator = ticker.MaxNLocator(nbins=7)
    #clbar.locator = tick_locator
    #clbar.update_ticks()
    clbar.set_label(colorbarLabelText) ###OFF for demonstration figure


    if False:  ### for demonstration figure
        y1 = 4996.0
        y2 = 5034.0
        y3 = 5090.0

        x1 = 6.0
        x2 = 7.0
        x3 = x1 + (x2-x1)/(y2-y1)*(y3-y1)
        x4 =-57
    
        plt.plot([x1, x3], [y1, y3], 'r--', linewidth=5)
        plt.plot([x1, x2], [y1, y2], 'k-', linewidth=5)
        plt.plot([x2, x4], [y2, y3], 'k-', label='3\sigma limits', linewidth=5)
    
        plt.plot([-125, 125], [y3, y3], 'xkcd:maroon', linestyle=':', label='3\sigma limits', linewidth=5)


    
    #plt.tight_layout()
    plt.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)
    #print(titleText+'.pdf')
    
    

    if createPDF:
        plt.savefig(titleText+'.pdf')
    else:
        plt.savefig(titleText+'.png')
        
    plt.show()
    
    #return(maximumValue)
    return(averageValue)

##########################################################################################################################################################
def sumHitForPanels(HitHist):

    DummyHistData, DummyHistAxisX, DummyHistAxisY, _, DummyLabelX, DummyLabelY = disassembleHist2d(HitHist)

    [nBinsXDummy,nBinsYDummy] = DummyHistData.shape

  
    nBinsYDummy
    
    panelSum=[0, 0, 0, 0, 0]
    for i in range(nBinsYDummy): #panel number
        for j in range(nBinsXDummy): #straw number

            # there was a change in storing data to histograms when moving to 40 tubes instead of 7
            # the tube number [0,40[ and straw number goes from  [0,7[ so the unique ID is (tube number)*100 + straw number
            # instead of the old (tube number)*10 + straw number
            if(nBinsYDummy<100): #
                panelNumber = int(i // 10)
            else:
                panelNumber = int(i // 100)
                
            if(panelNumber < 5):
                panelSum[panelNumber] += DummyHistData[j,i]
             
    print(panelSum)
    print(sum(panelSum))

    return panelSum



##########################################################################################################################################################
def plotHistogramFromList(valuesArray, yLabelText, yLimits=[0.0,1.0]):

    print(valuesArray)
    
    x = np.arange(5)
    print(x)
    fig, ax = plt.subplots()
    plt.bar(x, valuesArray)
    
    plt.xticks(x, ('1', '2', '3', '4' , '5'))
    #ax.set_ylim([0.28,0.35]) #lambda
    ax.set_ylim(yLimits)

    fig.set_size_inches(13.9,9.75)
    ax.tick_params(direction='in', bottom=1, top=1, left=1, right=1 )
    plt.tight_layout(pad=1.4)
    
    plt.xlabel('Panel')
    plt.ylabel(yLabelText)
    plt.title(' ')
    
    #plt.show()
    plt.savefig(yLabelText+'.pdf')
    
###########################################################################################################################################################

def decomposeFileName(fileName):
    parts=fileName.split('_')
    polyethyleneThickness = 0
    for part in parts:
        if 'ang' in part:
            wavelength=float(part.strip('ang').replace('p','.'))
        elif 'panel' in part:
            panelNumber = int(part.strip('panel'))
        elif 'with' in part:
            material = part.strip('with')
        elif 'PE' in part:
            polyethyleneThickness = int(part.strip('PE'))

    return [wavelength, panelNumber, material, polyethyleneThickness]
