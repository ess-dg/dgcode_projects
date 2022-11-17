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

# import BCS.BCSplotUtils as BCSplt
import LOKI.BCSplotUtils as BCSplt
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
def plotHardCodedConvLayerThicknessEfficiency(createPDF=1):

    
    values0p6 = [22.90, 27.29, 27.93]
    values1p8 = [45, 50, 49.6]
    values3 = [58.39, 59.64, 58.64]
    values5 = [67.81, 65.04, 63.04]
    values11 = [74.19, 66.30, 63.61]
    
    xAxisValues = [0.65,1.0,1.1]
    
    legendList = ['$\lambda$=0.6 {\AA}','$\lambda$=1.8 {\AA}', '$\lambda$=3 {\AA}', '$\lambda$=5 {\AA}', '$\lambda$=11 {\AA}']
    minValueList = list()
    maxValueList = list()

    fig = plt.figure()
    fig.set_size_inches(12,10)
    setPlottingParameters(1)
    #plt.rc('legend', fontsize=22)    # legend fontsize
    plt.rc('legend',**{'fontsize':30})
    plt.rc('legend',**{'labelspacing':0.3})
    
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for 5 wavelength
    
    mSize=10
    lWidth=4
    plt.plot(xAxisValues, values0p6, '--o' , linewidth=lWidth, markersize=mSize)
    plt.plot(xAxisValues, values1p8, '--o', linewidth=lWidth, markersize=mSize)
    plt.plot(xAxisValues, values3, '--o', linewidth=lWidth, markersize=mSize)
    plt.plot(xAxisValues, values5, '--o', linewidth=lWidth, markersize=mSize)
    plt.plot(xAxisValues, values11, '--o', linewidth=lWidth, markersize=mSize)
       
    plt.legend(legendList, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.03))
    plt.ylabel('Detection efficiency')
    plt.grid()

    yAxisMin = 0
    yAxisMax = 100
    xAxisMin = 0.5
    xAxisMax = 1.2
        
    plt.xlabel('Converter layer thickness [$\mu$m]')
    plt.axis([xAxisMin, xAxisMax, yAxisMin, yAxisMax])

    plt.plot([xAxisMin, xAxisMax], [70, 70], colors['magenta'], linewidth=4) #, label=limitsLabel)
    
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        plt.savefig('efficiencyConvThickness.pdf')
    else:
        plt.savefig('efficiencyConvThickness.png')
    plt.show()


 ##########################################################################################################################################################
def plotActivity(createPDF=1):
    
    setPlottingParameters(0.8)
    

    volumeAl = 14446.55 #cm^3
    volumeCu =  821.92 #cm^3
    volumeCu63 = volumeCu
    volumeCu65 = volumeCu
    
    #fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[1, 1]}, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[1,1]}, figsize=(16, 7), dpi=80, facecolor='w', edgecolor='k')


    ######## copper composition ##########
    abundanceCu63 = 0.6915
    abundanceCu65 = 0.3085
    sigmaTotalCu63 = 4.5
    sigmaTotalCu65 = 2.17
    absRatioCu63 = sigmaTotalCu63*abundanceCu63 / (sigmaTotalCu63*abundanceCu63 + sigmaTotalCu65*abundanceCu65)
    absRatioCu65 = 1 -absRatioCu63


    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
     
    ########### aluminium composition #########
    #used for aluminiumAlloy
    #sigmaAluminium = {
    #    'Al-28': [0.013114013, 135.00,  '#7f7f7f'],
    #    'Cr-51': [6.48066E-05, 2390000, colors['purple']],
    #    'Cr-55': [7.98668E-07, 210.00,  colors['indigo']],
    #    'Cu-64': [7.96166E-05, 45700.00, '#ff7f0e'],
    #    'Cu-66': [1.71283E-05, 307.00,   colors['orangered']],
    #    'Fe-55': [1.53178E-05, 86100000, '#DA1F28'],
    #    'Fe-59': [4.20425E-07, 3850000,  '#D7BB2B'],
    #    'Mg-27': [1.01312E-05, 567.00,   colors['darkred']],  
    #    'Mn-56': [0.001968102, 9280.00,   colors['magenta']],  
    #    'Si-31': [7.66031E-07, 9440.00,   colors['darkgreen']],
    #    'Ti-51': [4.72417E-07, 346.00,  colors['teal']],
    #    'Zn-65': [1.85843E-05, 21100000,  colors['royalblue']],
    #    'Zn-69': [9.1755E-05, 3380.00,    colors['cyan']],
    #    'Zn-71': [2.51792E-08, 147.00,    colors['navy']]}    
    #sigmaTotAl = 0.016403342 #This is more that the sum of the sigmaAluminium.values() because of the non active products

    sigmaAluminium = {'Al-28': [0.013114013, 135.00,  '#7f7f7f']} #used for pure aluminium
    sigmaTotAl = 0.013114013 #used for pure aluminium
    
    intensity = 1e9 * 0.05 #n/sec on the detector
    absRatioAl =  1.51/100 #3A
    absRatioCu =  1.76/100 #3A
    sourceAluminium = intensity * absRatioAl
    sourceCopper = intensity * absRatioCu
    
    sourceCu_63 = sourceCopper * absRatioCu63 #TODO
    sourceCu_65 = sourceCopper * absRatioCu65 #TODO

    sourceAl = {key:(sourceAluminium * sigmaAluminium[key][0]/sigmaTotAl) for key in sigmaAluminium.keys() }
    
    
    tIrrad = 1e6
    tCool = 2e6
   
    #x2 = np.linspace(tIrrad, tIrrad+tCool, 10000)
    
    #T12Al = 2.245 * 60 
    T12Cu_63 = 12.701 * 60*60 
    T12Cu_65 = 5.120 * 60 
    
    #lambdaAl = m.log(2)/T12Al
    lambdaCu_63 = m.log(2)/T12Cu_63
    lambdaCu_65 = m.log(2)/T12Cu_65
    lambdaAluminium = {key:(m.log(2) / sigmaAluminium[key][1]) for key in sigmaAluminium.keys() }
   

    #yEnd = 1e-10 #used for aluminiumAlloy
    yEnd = 1e-4 #used for pure aluminium

    yEND_Al = yEnd*volumeAl
    yEND_Cu63 = yEnd*volumeCu63
    yEND_Cu65 = yEnd*volumeCu65
    
    x = np.linspace(0.001, tIrrad, 20000)
    #xCool = np.linspace(1, tCool, 20000)
    
    
    #yAl = sourceAl * (1 - np.exp(-lambdaAl * x))
    yCu_63 = sourceCu_63 * (1 - np.exp(-lambdaCu_63 * x))
    yCu_65 = sourceCu_65 * (1 - np.exp(-lambdaCu_65 * x))
    yAluminium = {key:(sourceAl[key] * (1 - np.exp(-lambdaAluminium[key] * x))) for key in sigmaAluminium.keys() }
   
    #xEndAl = np.log(yEND_Al / yAl[-1]) / -lambdaAl
    xEndCu_63 = np.log(yEND_Cu63 / yCu_63[-1]) / -lambdaCu_63
    xEndCu_65 = np.log(yEND_Cu65 / yCu_65[-1]) / -lambdaCu_65
    xEndAl = {key: (np.log(yEND_Al / yAluminium[key][-1]) / -lambdaAluminium[key]) for key in sigmaAluminium.keys()}
    
    #xCoolAl = np.linspace(1, xEndAl, 20000)
    xCoolCu_63 = np.linspace(1, xEndCu_63, 20000)
    xCoolCu_65 = np.linspace(1, xEndCu_65, 20000)
    xCoolAluminium = {key: (np.linspace(1, xEndAl[key], 20000)) for key in sigmaAluminium.keys()}
    #xCoolCu_63 = xCoolAluminium['Cu-64'] #used for aluminiumAlloy
    #xCoolCu_65 = xCoolAluminium['Cu-66'] #used for aluminiumAlloy
    
    #yAlCool = yAl[-1] * np.exp(-lambdaAl * xCoolAl)
    yCu_64Cool = yCu_63[-1] * np.exp(-lambdaCu_63 * xCoolCu_63)
    yCu_65Cool = yCu_65[-1] * np.exp(-lambdaCu_65 * xCoolCu_65)
    yAluminiumCool = {key: (yAluminium[key][-1] * np.exp(-lambdaAluminium[key] * xCoolAluminium[key])) for key in sigmaAluminium.keys() }

    
    ############## Ploting ######################
    
    #ax1.plot(x, yAl/volumeAl, '#7f7f7f' , label = 'Al-27')
    ax1.plot(x, yCu_63/volumeCu63, '#ff7f0e', label = 'Cu-64')
    ax1.plot(x, yCu_65/volumeCu65, colors['orangered'], label = 'Cu-66')
    ax1.plot(x, yAluminium['Al-28']/volumeAl, '#7f7f7f' , label = 'Al-28')

    #ax2.plot(xCoolAl, yAlCool/volumeAl, '#7f7f7f', label = 'Al-28')
    ax2.plot(xCoolCu_63, yCu_64Cool/volumeCu63, '#ff7f0e', label = 'Cu-64')
    ax2.plot(xCoolCu_65, yCu_65Cool/volumeCu65, colors['orangered'], label = 'Cu-66')
    ax2.plot(xCoolAluminium['Al-28'], yAluminiumCool['Al-28']/volumeAl, '#7f7f7f', label = 'Al-28')

    decayGammaYield = {
        'Al-28':1.0 ,
        'Cu-64':0.353, #2.89e-3 , #1.22 with positrons
        'Cu-66':9.224e-2 }
    decayGammaIntensity = {
        'Al-28':decayGammaYield['Al-28']*yAluminium['Al-28'][-1]/volumeAl ,
        'Cu-64':decayGammaYield['Cu-64']*yCu_63[-1]/volumeCu63 ,
        'Cu-66':decayGammaYield['Cu-66']*yCu_65[-1]/volumeCu65 }

    print ('Cu volume: ' + str(volumeCu) + ' cm3');
    print ('Al volume: ' + str(volumeAl) + ' cm3');
    
    print ('Cu-64 saturation ativity concentration ' + str(yCu_63[-1]/volumeCu63)+' Bq')
    print ('Cu-66 saturation ativity concentration: ' + str(yCu_65[-1]/volumeCu65)+' Bq')
    print ('Al-28 saturation ativity concentration: ' + str(yAluminium['Al-28'][-1]/volumeAl )+' Bq')

    sumDecayGammaIntensity = (decayGammaIntensity['Cu-64']+decayGammaIntensity['Cu-66'])*volumeCu + decayGammaIntensity['Al-28']*volumeAl
    
    print ('Cu-64 decay gamma Intensity: ' + str(decayGammaIntensity['Cu-64'])+' 1/cm^3')
    print ('Cu-66 decay gamma Intensity: ' + str(decayGammaIntensity['Cu-66'])+' 1/cm^3')
    print ('Al-28 decay gamma Intensity: ' + str(decayGammaIntensity['Al-28'])+' 1/cm^3')
 
    print ('Cu decay gamma Intensity: ' + str((decayGammaIntensity['Cu-64']+decayGammaIntensity['Cu-66'])*volumeCu)+' 1/s')
    print ('Al decay gamma Intensity: ' + str(decayGammaIntensity['Al-28']*volumeAl)+' 1/s')
    print ('Sum decay gamma Intensity: ' + str(sumDecayGammaIntensity)+' 1/s')

    promptGammaYield = {
        'Al':1.978 ,
        'Cu':2.665 }
    promptGammaIntensity = {
        'Al':promptGammaYield['Al']*sourceAluminium ,
        'Cu':promptGammaYield['Cu']*sourceCopper}

    sumPromptGammaIntensity = promptGammaIntensity['Al']+promptGammaIntensity['Cu']

    print ('Cu prompt gamma intensity: ' + str(promptGammaIntensity['Cu']))
    print ('Al prompt gamma intensity: ' + str(promptGammaIntensity['Al']))
    print ('Sum prompt gamma intensity: ' + str(sumPromptGammaIntensity))

    print ('------------------------------')
    print ('Cu-64 decay gamma Intensity: ' + str((decayGammaIntensity['Cu-64'])*volumeCu)+' 1/s')
    print ('Cu-66 decay gamma Intensity: ' + str((decayGammaIntensity['Cu-66'])*volumeCu)+' 1/s')
    print ('Al decay gamma Intensity: ' + str(decayGammaIntensity['Al-28']*volumeAl)+' 1/s')
  
    ############# Editing figure #################
    
    #ax1.errorbar(alHistAxis, alHistData, yerr=alHistError, fmt='o-', label='Aluminium' )

    xLabel = 'time [s]'
    yLabel = 'Activity concentration [Bq/cm$^3$]'
    
    ax1.set_xlabel('Irradiation '+xLabel)
    ax1.set_ylabel(yLabel)
    ax1.legend(frameon=False)
    ax1.grid()
    ax1.set_ylim([1e-4, 2*1e3])
    ax1.set_xlim([0.01, tIrrad])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_xticklabels(["10$^{-2}$","10$^{-1}$","10$^0$","10$^1$", "10$^2$", "10$^3$", "10$^4$", "10$^5$", "10$^6$"])

    ax1Top = ax1.twiny()
    ax1Top.xaxis.set_ticks_position('top')
    ax1Top.set_xlim([0.01, tIrrad])
    ax1Top.set_xscale('log')
    ax1Top.set_xticks([60, 600, 60*60,24*60*60 ,7*24*60*60])
    ax1Top.set_xticklabels(["1 m","10 m","1 h","1 d","1 w"])
    
    
    ax2.set_xlabel('Cooling '+xLabel)
    #ax2.set_ylabel(yLabel)
    ax2.legend(frameon=False)
    ax2.grid()
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim([1e-4, 2*1e3])
    ax2.set_xlim([1, tCool])
    ax2.set_xticks([ 1e0, 1e1, 1e2, 1e3, 1e4,1e5,1e6])
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_xticklabels(["10$^0$","10$^1$", "10$^2$", "10$^3$", "10$^4$", "10$^5$", "10$^6$"])
   
    

    ax2Top = ax2.twiny()
    ax2Top.xaxis.set_ticks_position('top')
    ax2Top.set_xlim([1, tCool])
    ax2Top.set_xscale('log')

    ax2Top.set_xticks([60, 600, 60*60, 24*60*60,7*24*60*60])
    ax2Top.set_xticklabels(["1 m","10 m","1 h","1 d","1 w"])
    
    #plt.rc('legend',**{'fontsize':30})
    #plt.rc('legend',**{'labelspacing':0.3})
    
   
    #plt.legend(legendText, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.03))

    plt.subplots_adjust(wspace=0.05)
    fig.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)

    
    #ax1.set_yticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    #ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax1.set_yticklabels(["10$^0$","10$^1$", "10$^2$", "10$^3$", "10$^4$", "10$^5$", "10$^6$"])

    
    if createPDF:
        plt.savefig('activity.pdf')
    else:
        plt.savefig('activity.png')
    plt.show()


   ##########################################################################################################################################################
def plotPathLengthDistribution(fileName, createPDF=1):

    setPlottingParameters(0.8)
    xLabel = 'Path length [mm]'
    yLabel = 'Counts'
    
    hc = sh.HistCollection(fileName)
    alHist= hc.hist('neutron_segment_length_TubeWall')
    cuHist= hc.hist('neutron_segment_length_StrawWall')
    
    [alHistData, alHistError, alHistAxis, alTitleText, alHistLabel] = disassembleHist1d(alHist)
    [cuHistData, cuHistError, cuHistAxis, cuTitleText, cuHistLabel] = disassembleHist1d(cuHist)

    #fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[1, 1]}, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    fig, (ax1, ax2) = plt.subplots(2,1, gridspec_kw = {'width_ratios':[1]}, figsize=(14, 8), dpi=80, facecolor='w', edgecolor='k')

    ax1.errorbar(alHistAxis, alHistData, yerr=alHistError, fmt='o-', label='Aluminium' )
    ax2.errorbar(cuHistAxis, cuHistData, yerr=cuHistError, fmt='o-', label='Copper' )

    #ax1.set_xlabel(xLabel)
    ax2.set_xlabel(xLabel)

    ax1.set_ylabel(yLabel)
    ax2.set_ylabel(yLabel)

    
    ax1.legend(frameon=False)
    ax2.legend(frameon=False)
    
    ax1.grid()
    ax2.grid()
    #ax2.yaxis.tick_right()
    
    #ax2.set_ylim([yM, yMax2]) 
    ax1.set_xlim([0, 100])
    ax2.set_xlim([0, 5])

    ax1.set_yscale('log')
    ax2.set_yscale('log')
    

    
    #plt.rc('legend',**{'fontsize':30})
    #plt.rc('legend',**{'labelspacing':0.3})
    
   
    #plt.legend(legendText, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.03))

    plt.subplots_adjust(wspace=0.15)
    fig.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)

    
    ax1.set_yticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_yticklabels(["10$^0$","10$^1$", "10$^2$", "10$^3$", "10$^4$", "10$^5$", "10$^6$"])

    ax2.set_yticks([1e0,1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_yticklabels(["10$^0$","10$^1$", "10$^2$", "10$^3$", "10$^4$", "10$^5$", "10$^6$"])

    if createPDF:
        plt.savefig('pathLength.pdf')
    else:
        plt.savefig('pathLength.png')
    plt.show()


##########################################################################################################################################################
def plotHardCodedChart(createPDF=1):

   
    legends = ('1', '2', '3', '4', '5')
    legends = ('Panel 1', 'Panel 2', 'Panel 3', 'Panel 4', 'Panel 5')
    ticks=['0.6','1.8', '3', '5', '11']     
        
    meanMatrix = [[23.84, 22.38, 20.18, 17.94, 15.66],
                  [31.59, 25.09, 18.75, 13.98, 10.59],
                  [38.51, 26.95, 16.86, 10.66, 7.01],
                  [48.21, 28.18, 13.42, 6.66, 3.53],
                  [68.66, 24.87, 4.92, 1.24, 0.33]]

    nHit0p6 = 5.45586e6
    nHit1p8 = 9.91854e6
    nHit3 = 1.19182e7
    nHit5 = 1.29957e7
    nHit11 = 1.3270907e7

    def calcError(x,n):
        return m.sqrt((x**2 + 100*x)/(1e4*n))
    
    stdMatrix = [[calcError(23.84,nHit0p6), calcError(22.38,nHit0p6), calcError(20.18,nHit0p6), calcError(17.94,nHit0p6), calcError(15.66,nHit0p6)],
                  [calcError(31.59,nHit1p8), calcError(25.09,nHit1p8), calcError(18.75,nHit1p8), calcError(13.98,nHit1p8), calcError(10.59,nHit1p8)],
                  [calcError(38.51,nHit3), calcError(26.95,nHit3), calcError(16.86,nHit3), calcError(10.66,nHit3), calcError(7.01,nHit3)],
                  [calcError(48.21,nHit5), calcError(28.18,nHit5), calcError(13.42,nHit5), calcError(6.66,nHit5), calcError(3.53,nHit5)],
                  [calcError(68.66,nHit11), calcError(24.87,nHit11), calcError(4.92,nHit11), calcError(1.24,nHit11), calcError(0.33,nHit11)]]

                                                        

    plotBarChart(meanMatrix, stdMatrix, ticks, legends, createPDF=createPDF)

##########################################################################################################################################################
def plotHardCodedPanelEfficiency(createPDF=1):

    xAxisQantity = 'panelNumber'
    
    values0p6 = [6.50, 12.61, 18.11, 23.00, 27.27]
    values1p8 = [15.24, 27.74, 36.98, 44.17, 49.56]
    values3 = [22.95, 39.01, 49.06, 55.41, 59.59]
    values5 = [31.36, 49.69, 58.41, 62.75, 65.04]
    values11 = [45.57, 62.07, 65.34, 66.16, 66.37]

    def calcError(x):
        N=2e7
        return m.sqrt(x*100/N) 
    
    errors0p6 = [calcError(6.50), calcError(12.61), calcError(18.11), calcError(23.00), calcError(27.27)]
    errors1p8 = [calcError(15.24), calcError(27.74), calcError(36.98), calcError(44.17), calcError(49.56)]
    errors3 = [calcError(22.95), calcError(39.01), calcError(49.06), calcError(55.41), calcError(59.59)]
    errors5 = [calcError(31.36), calcError(49.69), calcError(58.41), calcError(62.75), calcError(65.04)]
    errors11 = [calcError(45.57), calcError(62.07), calcError(65.34), calcError(66.16), calcError(66.37)]
    
    
    xAxisValues = [1,2,3,4,5]
    
    legendList = ['$\lambda$=0.6 {\AA}','$\lambda$=1.8 {\AA}', '$\lambda$=3 {\AA}', '$\lambda$=5 {\AA}', '$\lambda$=11 {\AA}']
    minValueList = list()
    maxValueList = list()

    fig = plt.figure()
    #fig.set_size_inches(12,10.2)
    #setPlottingParameters(1)
    #plt.rc('legend',**{'fontsize':30})
    #plt.rc('legend',**{'labelspacing':0.3})
    fig.set_size_inches(8, 6)
    setPlottingParameters(0.8)
    plt.tight_layout()
    
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for 5 wavelength
    
    mSize=10
    lWidth=4

    plt.errorbar(xAxisValues, values0p6, yerr=errors0p6, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values1p8, yerr=errors1p8, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values3, yerr=errors3, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values5, yerr=errors5, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values11, yerr=errors11, fmt='o--', linewidth=lWidth, markersize=mSize)
       
    plt.legend(legendList, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.03))
    plt.ylabel('Global detection efficiency [\%]')
    plt.grid()

    yAxisMin = 0
    yAxisMax = 30
    xAxisMin = 0.0
    xAxisMax = 0.4
        
    plt.xlabel('Panel number')
    plt.axis([xAxisMin, xAxisMax, yAxisMin, yAxisMax])

    plt.plot([xAxisMin, xAxisMax], [70, 70], colors['magenta'], linewidth=4) #, label=limitsLabel)
    
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        plt.savefig('efficiencyPanel.pdf')
    else:
        plt.savefig('efficiencyPanel.png')
    plt.show()

 
    
##########################################################################################################################################################
def plotSignalChart(fileName, counterNames=0,createPDF=0):

    legends = (['Signal'])
    if not counterNames:
        counterNames = ['dx', 'dy', 'dth', 'dphi', 'dtof', 'dlambda', 'dq']

    counterLabel ={
        "dx": "$\delta$X",
        "dy": "$\delta$Y",
        "dth": "$\delta\Theta$",
        "dphi": "$\delta\Phi$",
        "dtof": "$\delta$TOF",
        "dlambda": "$\delta\lambda$",
        "dq": "$\delta$Q",
    }
        
    meanMatrix = list()
    stdMatrix = list()

    ticks=list()     

    signalMin=100.0
    signalMax=0.0
    
    hc = sh.HistCollection(fileName)
    for counter in counterNames:
        numberOfNeutrons = float(hc.hist('neutron_true_x').getIntegral())
        
        signal = float(hc.hist("signal_counters").getCounter(counter).getValue())/ numberOfNeutrons
        signalError = float(hc.hist("signal_counters").getCounter(counter).getError())/ numberOfNeutrons
        
        meanMatrix.append([signal*100])
        signalMin = min(signalMin, signal*100)
        signalMax = max(signalMax, signal*100)
        
        stdMatrix.append([signalError*100])
        
        ticks.append(counterLabel[counter])

    print (str( '%.2f' % (signalMax - signalMin)) + ' File: '+fileName)
    print ('min: ' + str(signalMin))
    print ('max: ' + str(signalMax))
    setPlottingParameters(0.6)
    plotBarChart(meanMatrix, stdMatrix, ticks, legends, createPDF=createPDF, xLabelText='Quantity')

##########################################################################################################################################################
def plotAbsorptionChart(fileNames, counterName='dth', plotSignal=0 ,createPDF=0): #'dth'

    if plotSignal:
        legends = ('Signal', 'Background', 'Transmission', 'Absorption')
    else:
        legends = ('Detection', 'Transmission', 'Abs in B4C', 'Abs in Cu', 'Abs in Al')

    meanMatrix = list()
    stdMatrix = list()

    ticks=list()     

    for i in range(len(fileNames)):
        [wavelength, panelNumber, material, polyethyleneThickness] = decomposeFileName(fileNames[i])
        
        hc = sh.HistCollection(fileNames[i])
        numberOfNeutrons = float(hc.hist('neutron_true_x').getIntegral())

        signal = float(hc.hist("signal_counters").getCounter(counterName).getValue())/ numberOfNeutrons
        signalError = float(hc.hist("signal_counters").getCounter(counterName).getError())/ numberOfNeutrons
        
        detection = float((hc.hist('neutron_x_hit')).getIntegral()) / numberOfNeutrons
        detectionError = m.sqrt(float((hc.hist('neutron_x_hit')).getIntegral())) / numberOfNeutrons
        background = detection - signal
        backgroundError = m.sqrt(detectionError**2 + signalError**2)
        
        converter = float((hc.hist('neutron_x_conv')).getIntegral()) / numberOfNeutrons #neutrons that end in B4C
        converterError = m.sqrt(float((hc.hist('neutron_x_conv')).getIntegral()))/ numberOfNeutrons
        boron_abs = converter - detection
        boron_absError = m.sqrt(converterError**2 + detectionError**2)
        tubeWall_hit = float(hc.hist('h_st_tubeWall_hit').getIntegral()) / numberOfNeutrons
        tubeWall_hitError = m.sqrt(float(hc.hist('h_st_tubeWall_hit').getIntegral())) / numberOfNeutrons
        strawWall_hit = float(hc.hist('h_st_strawWall_hit').getIntegral()) / numberOfNeutrons
        strawWall_hitError = m.sqrt(float(hc.hist('h_st_strawWall_hit').getIntegral())) / numberOfNeutrons
        absorption = boron_abs + tubeWall_hit + strawWall_hit
        absorptionError = m.sqrt(boron_absError**2 + strawWall_hitError**2)
        
        transmission = 1.0 - (detection + absorption)
        transmissionError = m.sqrt(detectionError**2 + absorptionError**2)
        
        if plotSignal:
            meanMatrix.append([signal*100, background*100, transmission*100, absorption*100])
            stdMatrix.append([signalError*100, backgroundError*100, transmissionError*100, absorptionError*100])

            print([signal*100, background*100, transmission*100, absorption*100])
            
            if polyethyleneThickness :
                tictTextPE = '\nPE'
                emptyColumnForSeparation = 1
            else:
                tictTextPE = ''
                emptyColumnForSeparation = 0
            ticks.append(str(wavelength)+tictTextPE)
            

            if emptyColumnForSeparation and i != len(fileNames)-1:
                meanMatrix.append([0, 0, 0, 0])
                stdMatrix.append([0, 0, 0, 0])
                ticks.append('')
        else:
            meanMatrix.append([detection*100, transmission*100, boron_abs*100, strawWall_hit*100, tubeWall_hit*100])
            stdMatrix.append([detectionError*100, transmissionError*100, boron_absError*100, strawWall_hitError*100, tubeWall_hitError*100])

            #print [detection*100, transmission*100, boron_abs*100, strawWall_hit*100, tubeWall_hit*100]
            
            if polyethyleneThickness :
                tictTextPE = '\nPE'
                emptyColumnForSeparation = 1
            else:
                tictTextPE = ''
                emptyColumnForSeparation = 0
            ticks.append(str(wavelength)+tictTextPE)

            if emptyColumnForSeparation and i != len(fileNames)-1:
                meanMatrix.append([0, 0, 0, 0, 0])
                stdMatrix.append([0, 0, 0, 0, 0])
                ticks.append('')
        

    plotBarChart(meanMatrix, stdMatrix, ticks, legends, createPDF=createPDF)
    
##########################################################################################################################################################
def plotBarChart(meanMatrix, stdMatrix, ticks, legends, createPDF=0, xLabelText='Wavelength [{\AA}]'):

    N = len(ticks)
    M =  len(legends)
    
    mean = [ np.array(meanMatrix)[:,i] for i in reversed(list(range(M)))]
    std = [ np.array(stdMatrix)[:,i] for i in reversed(list(range(M)))]

    buttom = list()  
    buttom.append(tuple([0 for i in range(N)]))
    for i in range(M-1):
        buttom.append(tuple(map(operator.add, buttom[i], mean[-(1+i)])))
    
    ind = np.arange(N)    # the x locations for the groups
    width =0.6 # 0.9     # the width of the bars: can also be len(x) sequence
    #width =0.9 # PE

    fig = plt.figure()
    fig.set_size_inches(8, 6)
    setPlottingParameters(0.8)
    plt.tight_layout()

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    if(len(legends)==4 or len(legends)==1):
        plt.rc('axes', prop_cycle=(cycler('color', [ '#1f77b4',  '#d62728', 'gold','rebeccapurple'])))  #green '#2ca02c'
    elif(len(legends)==5):
        plt.rc('axes', prop_cycle=(cycler('color', [ 'green', 'gold' , (0.0, 0.0, 0.0, 0.89), '#ff7f0e', '#7f7f7f']))) #blue '#1f77b4'
    else:
      print('todo set color')

    p = [ plt.bar(ind, mean[-(1+i)], width, bottom=buttom[i], yerr=std[-(1+i)])
        for i in range(M)]
    

    #plt.ylabel('Relative ratio [{\%}]')
    plt.ylabel('Relative fraction [{\%}]')
    #plt.ylabel('Relative signal fraction [{\%}]')
    #plt.title('Wavelength')
    plt.xlabel(xLabelText)
    plt.xticks(ind, ticks)
    #TEMPplt.yticks(np.arange(0, 110, 10))

    #plt.plot([-1, 7], [44.06597, 44.06597], 'r:', linewidth=4)
    #plt.plot([-1, 7], [46.910035, 46.910035], 'r:', linewidth=4)
    #plt.xlim(-0.7,6.7) #noAvoid = False
    plt.xlim(-0.7,3.7) #noAvoid = True
    
    plt.legend([p[i][0] for i in range(M)], legends, frameon=True, loc='lower right') # loc='top right') # 

    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        plt.savefig('barChart'+legends[0]+'.pdf')
    else:
        plt.savefig('barChart'+legends[0]+'.png')
    plt.show()

##########################################################################################################################################################
def fractionalScatteringPlot(fileNames, createPDF=0):

    counterName = 'dth'
    histName = 'neutron_dth_hit'
    xAxisQantity = 'wavelength'

    valuesAl = list()
    valuesCu = list()
    valuesAlCu = list()
    valuesB = list()
    valuesPE = list()
    errorsAl = list()
    errorsCu = list()
    errorsAlCu = list()
    errorsB = list()
    errorsPE = list()
    xAxisValuesAl = list()
    xAxisValuesCu = list()
    xAxisValuesAlCu = list()
    xAxisValuesB = list()
    xAxisValuesPE = list()
    legendList = list()
    minValueList = list()
    maxValueList = list()

    def calcError(x,n):
        return m.sqrt((x**2 + x*n)/(n**3))
   
    for i in range(len(fileNames)):
        [wavelength, panelNumber, material, polyethyleneThickness] = decomposeFileName(fileNames[i])
        
        hc = sh.HistCollection(fileNames[i])
        #print hc.hist("signal_counters").getCounter(counterName).getValue()
        signal = float(hc.hist("signal_counters").getCounter(counterName).getValue())
        signalAndBackground = float((hc.hist(histName)).getIntegral())
        #sumAbsorption = tubeWall_hit + strawWall_hit
        ##numberOfNeutrons = float(hc.hist('neutron_true_x').getIntegral())
        #numberOfHits = float(hc.hist('neutron_x_hit').getIntegral())
        fractionalScattering = (signalAndBackground - signal) / signalAndBackground
        

        if material == 'Al':
            valuesAl.append(fractionalScattering)
            errorsAl.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValuesAl.append(wavelength)
            
        elif material == 'Cu':
            valuesCu.append(fractionalScattering)
            errorsCu.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValuesCu.append(wavelength)
             
        elif material == 'AlCu':
            if polyethyleneThickness : #simulation with PE
                valuesPE.append(fractionalScattering)
                errorsPE.append(calcError((signalAndBackground - signal), signalAndBackground))
                xAxisValuesPE.append(wavelength)
            else:
                valuesAlCu.append(fractionalScattering)
                errorsAlCu.append(calcError((signalAndBackground - signal), signalAndBackground))
                xAxisValuesAlCu.append(wavelength)
                
          
        elif material == 'B':
            valuesB.append(fractionalScattering)
            errorsB.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValuesB.append(wavelength)
        
        else:
            print('Problem with the material name in the filename of: '+ fileNames[i])
            
    #fig = plt.figure()
    #ax = plt.gca()
    setPlottingParameters(0.6)

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']])))

    #plt.errorbar(xAxisValues0p6, values0p6, yerr=errors0p6, fmt='o--', linewidth=2, markersize=6)
                           
    maxAl = maxCu = maxAlCu = 0
    minAl = minCu = minAlCu = 0
    if xAxisValuesPE:
        plt.errorbar(xAxisValuesPE, valuesPE, yerr=errorsPE, fmt='o--', linewidth=2, markersize=6)
        legendList.append('B+Al+Cu+PE')
        maxValueList.append(max(valuesPE))
        minValueList.append(min(valuesPE))
    if xAxisValuesAlCu:
        plt.errorbar(xAxisValuesAlCu, valuesAlCu, yerr=errorsAlCu, fmt='o--', linewidth=2, markersize=6)
        legendList.append('B+Al+Cu')
        maxValueList.append(max(valuesAlCu))
        minValueList.append(min(valuesAlCu))
    if xAxisValuesAl:
        plt.errorbar(xAxisValuesAl, valuesAl, yerr=errorsAl, fmt='o--', linewidth=2, markersize=6)
        legendList.append('B+Al')
        maxValueList.append(max(valuesAl))
        minValueList.append(min(valuesAl))
    if xAxisValuesCu:
        plt.errorbar(xAxisValuesCu, valuesCu, yerr=errorsCu, fmt='o--', linewidth=2, markersize=6)
        legendList.append('B+Cu')
        maxValueList.append(max(valuesCu))
        minValueList.append(min(valuesCu))
    if xAxisValuesB:
        plt.errorbar(xAxisValuesB, valuesB, yerr=errorsB, fmt='o--', linewidth=2, markersize=6)
        legendList.append('B')
        maxValueList.append(max(valuesB))
        minValueList.append(min(valuesB))
    

    #plt.title('Fractional Scattering ('+re.escape(histName)+')')
    plt.legend(legendList, frameon=False)
    plt.ylabel('Fractional scattering')
    plt.grid()

    yAxisEmptySpace = 0.1*(max(maxValueList) - min(minValueList))
    yAxisMin = min(minValueList) - yAxisEmptySpace
    yAxisMax = max(maxValueList) + yAxisEmptySpace
   
    plt.xlabel('Wavelength [{\AA}]')
    plt.axis([0, 12, yAxisMin, yAxisMax])
    
    #fig.set_size_inches(13.9/1.8,9.75/1.8) #theGoodOne
    
    
    #fig.set_size_inches(2*6.3889,2*4.7917)
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        #plt.savefig('fracScatDTheta_MatWave.pdf')
        plt.savefig('fracScatDTheta_ang1p8_Panel.pdf')
    else:
        #plt.savefig('fracScatDTheta_MatWave.png')
        plt.savefig('fracScatDTheta_ang1p8_Panel.png')
    plt.show()
    
    ##########################################################################################################################################################
def fractionalScatteringPanelNrPlot(fileNames, createPDF=0):

    counterName = 'dth'
    histName = 'neutron_dth_hit'
    #xAxisQantity = 'wavelength'
    xAxisQantity = 'panelNumber'

    values0p6 = list()
    values1p8 = list()
    #values3 = list()
    values11 = list()
    errors0p6 = list()
    errors1p8 = list()
    errors11 = list()
    xAxisValues0p6 = list()
    xAxisValues1p8 = list()
    #xAxisValues3 = list()
    xAxisValues11 = list()
    legendList = list()
    minValueList = list()
    maxValueList = list()

    def calcError(x,n):
        return m.sqrt((x**2 + x*n)/(n**3))
    
    for i in range(len(fileNames)):
        [wavelength, panelNumber, material, ommit] = decomposeFileName(fileNames[i])
        
        hc = sh.HistCollection(fileNames[i])
        #print hc.hist("signal_counters").getCounter(counterName).getValue()
        signal = float(hc.hist("signal_counters").getCounter(counterName).getValue())
        signalAndBackground = float((hc.hist(histName)).getIntegral())
        #sumAbsorption = tubeWall_hit + strawWall_hit
        ##numberOfNeutrons = float(hc.hist('neutron_true_x').getIntegral())
        #numberOfHits = float(hc.hist('neutron_x_hit').getIntegral())
        fractionalScattering = (signalAndBackground - signal) / signalAndBackground
        

        if wavelength == 0.6:
            values0p6.append(fractionalScattering)
            errors0p6.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValues0p6.append(panelNumber)
        elif wavelength == 1.8:
            values1p8.append(fractionalScattering)
            errors1p8.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValues1p8.append(panelNumber)
        #elif wavelength == 3:
        #    values3.append(fractionalScattering)
        #    xAxisValues3.append(panelNumber)
        elif wavelength == 11:
            values11.append(fractionalScattering)
            errors11.append(calcError((signalAndBackground - signal), signalAndBackground))
            xAxisValues11.append(panelNumber)
        else:
            print('Problem with the wavelength in the filename of: '+ fileNames[i])
            
    #fig = plt.figure()
    #ax = plt.gca()
    setPlottingParameters(0.6)
    
    
    max0p6 = max1p8 = max3 = max11 = 0
    min0p6 = min1p8 = min3 = min11 = 0

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['navy']]))) #colors['darkgreen'], colors['cyan'],


    #plt.errorbar(wavelengthB, valuesB, yerr=errorB, fmt='o--', linewidth=2, markersize=6)
    
    if xAxisValues0p6:
        plt.errorbar(xAxisValues0p6, values0p6, yerr=errors0p6, fmt='o--', linewidth=2, markersize=6)
        legendList.append('0.6 {\AA}')
        maxValueList.append(max(values0p6))
        minValueList.append(min(values0p6))
    if xAxisValues1p8:
        plt.errorbar(xAxisValues1p8, values1p8, yerr=errors1p8, fmt='o--',linewidth=2, markersize=6)
        legendList.append('1.8 {\AA}')
        maxValueList.append(max(values1p8))
        minValueList.append(min(values1p8))
    #if xAxisValues3:
    #    plt.plot(xAxisValues3, values3, '--o', linewidth=2, markersize=6)
    #    legendList.append('3.0 [{\AA}]')
    #    maxValueList.append(max(values3))
    #    minValueList.append(min(values3))
    if xAxisValues11:
        plt.errorbar(xAxisValues11, values11, yerr=errors11, fmt='o--', linewidth=2, markersize=6)
        legendList.append('11.0 {\AA}')
        maxValueList.append(max(values11))
        minValueList.append(min(values11))

    #plt.title('Fractional Scattering ('+re.escape(histName)+')')
    plt.legend(legendList, frameon=False)
    plt.ylabel('Fractional scattering')
    plt.grid()

    yAxisEmptySpace = 0.1*(max(maxValueList) - min(minValueList))
    yAxisMin = min(minValueList) - yAxisEmptySpace
    yAxisMax = max(maxValueList) + yAxisEmptySpace
   

    plt.xlabel('Panel number')
    plt.axis([0.5, 5.5, yAxisMin, yAxisMax])
    #ax.set_ylim([yAxisMin, yAxisMax ])

    #fig.set_size_inches(13.9/1.8,9.75/1.8) #theGoodOne
    
    
    #fig.set_size_inches(2*6.3889,2*4.7917)
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        #plt.savefig('fracScatDTheta_MatWave.pdf')
        plt.savefig('fracScatDTheta_Panel.pdf')
    else:
        #plt.savefig('fracScatDTheta_MatWave.png')
        plt.savefig('fracScatDTheta_Panel.png')
    plt.show()
    
    
    ##########################################################################################################################################################
def fitGaussian(shistFileName, histName, createPDF=0, normalise=1, plotFit=1, plotFitAll=0, plotLimits=1, plotCenter=0,):
    """The gaussian fitting might fail because of bad initial mean and sigma values (0,0.1 worked in almost all cases this function has been used so far) """

    setPlottingParameters(0.8)
    
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    if 1 < len(shistFileName):
        plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for 5 wavelength
        #plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for material plots
        #plt.rc('axes', prop_cycle=(cycler('color', [colors['blue'],  colors['darkred']]))) #for PE plot 0.6Ang
    else:
        plt.rc('axes', prop_cycle=(cycler('color', [colors['black'],  colors['blue'], colors['darkred']])))
    plt.enable_tex_fonts()
    
    if plotCenter:
        fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[2, 1]}, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    else: 
        fig, ax1 = plt.subplots(1,1 , figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')

    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2] = getHardcodedLimits(histName[0], 11.0, 0, 0) #just to get xMin2, xMax2 !!
    xMinList = list()
    xMaxList = list()
    yMinList = list()
    yMaxList = list()
    yMinList2 = list()
    yMaxList2 = list()
    legendList = list()

    tempVar = 0

    for fileIndex in range(len(shistFileName)):
        for histIndex in range(len(histName)):
            hc = sh.HistCollection(shistFileName[fileIndex])
            my_histogram = hc.hist(histName[histIndex])

            #if(fileIndex==0): #Rebinnig needed for PE results
            #if(fileIndex>0): #Rebinnig needed  compare panel/material
            #    my_histogram.rebin(500)
            
            [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)

           
            x = np.array(histAxis)
            yErr = np.array(histError)
            normFactor = 1
            if normalise:
                normFactor = sum(histData)*(x[2]-x[1])
            y = [ tmp/normFactor for tmp in histData ]

            
            [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])

            if len(histName)>1:
                labelText = re.escape(str(histName[histIndex]))
            elif False:
                if material == 'AlCu':
                    labelText = 'B+Al+Cu'
                elif material == 'Al':
                    labelText = 'B+Al'
                elif material == 'Cu':
                    labelText = 'B+Cu'
                elif material == 'B':
                    labelText = 'B'
            elif False:
                if tempVar == 0:
                    tempVar = 1
                    labelText = 'B+Al+Cu+PE'
                else:
                    labelText = 'B+Al+Cu'
            elif False:
                labelText = 'Number of panels = '+str(panelNumber)
            else:
                labelText = '$\lambda$='+str(wavelength)+' {\AA}'
                
            ax1.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText )

            yNPArray = np.array(y)
            xNPArray = np.array(x)
            
            yMinList.append(min(yNPArray[np.nonzero(yNPArray)])) 
            yMaxList.append(max(y))
            
            xMinList.append(min(xNPArray[np.nonzero(yNPArray)]))
            xMaxList.append(max(xNPArray[np.nonzero(yNPArray)]))
            
            if plotCenter:
                ax2.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText, markersize=10 )
                tempYCenterList = list()
                for i in range(len(yNPArray)): #TODO CLEAN IT!
                    if xMin2 <= xNPArray[i] and xNPArray[i] <= xMax2 and yNPArray[i]>0.0 : 
                        tempYCenterList.append(yNPArray[i])
                yMinList2.append(min(tempYCenterList))
                yMaxList2.append(max(tempYCenterList))

        
            if plotFit==1 or plotFitAll:
                #estimate mean and standard deviation
                mean = 0 # np.dot(x,y)
                sigma = 0.1 #this works for almost all 
                #sigma = 0.001# USE this if the fitting fails
                #do the fit!
                popt, pcov = curve_fit(gauss_function, x, y, p0 = [1, mean, sigma])
                ### plot fitted curve ###
                xFit = x
                #plt.semilogy(x, gauss_function(x, *popt), 'b-', label='fitted curve')
                meanFit = popt[1]
                sigmaFit = popt[2]

                print(shistFileName[fileIndex], histName[histIndex],'Mean: ', popt[1] ,  'Sigma ' ,  popt[2])
                plotFit=42 #Plot fitted curve only for the first
       
    dx=x[1]-x[0]

    yMin = min(yMinList)/2
    yMax = max(yMaxList)*2
    xMin = min(xMinList)- 10 * dx
    xMax = max(xMaxList)+ 10 * dx

    #yMin2, Ymax2, xMin2, xMax2 = yMin, yMax, xMin, xMax
    
    limitsLabel = 'Signal limits'
    
    
    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2 ] = getHardcodedLimits(histName[0], wavelength, xMin, xMax)        
    ax1.set_xlabel(xLabel)
    #fig.text(0.5, 0.04, xLabel, ha='center')

    ax1.set_ylabel('Counts')
    
    
    if plotFit!=0:
        ax1.semilogy(x, gauss_function(x, *popt), '-', label='Fitted curve', zorder=10)
         
    if plotLimits:
        if plotFit!=0 and gaussian==1: #else the limits are already set earlier
            limitLeft = meanFit - 3 * sigmaFit
            limitRight = meanFit + 3 * sigmaFit
        
        ax1.plot([limitLeft, limitLeft], [yMin, yMax], colors['magenta'])
        ax1.plot([limitRight, limitRight], [yMin, yMax], colors['magenta'], label=limitsLabel)
 

    ax1.set_ylim([yMin, yMax]) 
    ax1.set_xlim([xMin, xMax])
    ax1.set_yscale('log')
    #ax1.legend(frameon=False, loc='upper left', bbox_to_anchor=(-0.01, 1.03))
    ax1.legend(frameon=False)
    ax1.grid()

    if plotCenter:
        yMin2 = min(yMinList2)/2
        yMax2 = max(yMaxList2)*2 
        if plotFit!=0:
            ax2.semilogy(x, gauss_function(x, *popt), 'b-', label='Fitted curve', zorder=10)
        if plotLimits:
            ax2.plot([limitLeft, limitLeft], [yMin2, yMax2], colors['magenta'])
            ax2.plot([limitRight, limitRight], [yMin2, yMax2], colors['magenta'], label=limitsLabel)
            
        ax2.set_xlabel(xLabel)
        #ax2.yaxis.tick_right() #is it better?
        
        ax2.set_ylim([yMin2, yMax2]) 
        ax2.set_xlim([xMin2, xMax2])
        ax2.set_yscale('log')
        ax2.grid()
    
    fig.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)
    #fig.tight_layout(pad=0.4, w_pad=0.2, h_pad=0.1)
    
    if createPDF:
        fig.savefig('gaussianFit_'+histName[0]+'.pdf')
    else:
        fig.savefig('gaussianFit_'+histName[0]+'.png')
    fig.show()


    ##########################################################################################################################################################
def fitGaussianEKIN(shistFileName, histName, createPDF=0, normalise=1, plotFit=1, plotFitAll=0, plotLimits=1, plotCenter=0,):
    """The gaussian fitting might fail because of bad initial mean and sigma values (0,0.1 worked in almost all cases this function has been used so far) """

    setPlottingParameters(0.7)
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    if 1 < len(shistFileName):
        plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for 5 wavelength
        #plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for material plots
        #plt.rc('axes', prop_cycle=(cycler('color', [colors['blue'],  colors['darkred']]))) #for PE plot 0.6Ang
    else:
        plt.rc('axes', prop_cycle=(cycler('color', [colors['black'],  colors['blue'], colors['darkred']])))
    plt.enable_tex_fonts()
    
    if plotCenter:
        fig, (ax1, ax2, ax3) = plt.subplots(1,3, gridspec_kw = {'width_ratios':[2, 1, 1]}, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    else: 
        fig, ax1 = plt.subplots(1,1 , figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')

    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2] = getHardcodedLimits(histName[0], 11.0, 0, 0) #just to get xMin2, xMax2 !!
    xMinList = list()
    xMaxList = list()
    yMinList = list()
    yMaxList = list()
    yMinList2 = list()
    yMaxList2 = list()
    yMinList3 = list()
    yMaxList3 = list()
    legendList = list()

    tempVar = 0

    for fileIndex in range(len(shistFileName)):
        for histIndex in range(len(histName)):
            print(shistFileName[fileIndex])
            hc = sh.HistCollection(shistFileName[fileIndex])
            my_histogram = hc.hist(histName[histIndex])
    
            [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)

            x = np.array(histAxis)
            yErr = np.array(histError)
            normFactor = 1
            if normalise:
                normFactor = sum(histData)*(x[2]-x[1])
            y = [ tmp/normFactor for tmp in histData ]

            
            [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])


            labelText = '$\lambda$='+str(wavelength)+' {\AA}'
                
            ax1.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText )

            yNPArray = np.array(y)
            xNPArray = np.array(x)
            
            yMinList.append(min(yNPArray[np.nonzero(yNPArray)])) 
            yMaxList.append(max(y))
            
            xMinList.append(min(xNPArray[np.nonzero(yNPArray)]))
            xMaxList.append(max(xNPArray[np.nonzero(yNPArray)]))
            
            if plotCenter:

                my_histogram = hc.hist(histName[histIndex]+'_center')
    
                [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)
            
                x = np.array(histAxis)
                yErr = np.array(histError)
                normFactor = 1
                if normalise:
                    normFactor = sum(histData)*(x[2]-x[1])
                y = [ tmp/normFactor for tmp in histData ]

                [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])

                yNPArray = np.array(y)
                xNPArray = np.array(x)
        
                ax2.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText , markersize=10 )
                tempYCenterList = list()
                for i in range(len(yNPArray)): #TODO CLEAN IT!
                    if xMin2 <= xNPArray[i] and xNPArray[i] <= xMax2 and yNPArray[i]>0.0 : 
                        tempYCenterList.append(yNPArray[i])
                yMinList2.append(min(tempYCenterList))
                yMaxList2.append(max(tempYCenterList))
                
                #################################################################################
                my_histogram = hc.hist(histName[histIndex]+'_center2')
    
                [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)
            
                x = np.array(histAxis)
                yErr = np.array(histError)
                normFactor = 1
                if normalise:
                    normFactor = sum(histData)*(x[2]-x[1])
                y = [ tmp/normFactor for tmp in histData ]

                [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])

                yNPArray = np.array(y)
                xNPArray = np.array(x)
                
                xMin3, xMax3 = -0.05, 0.05                
                ax3.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText , markersize=10 )
                tempYCenterList = list()
                for i in range(len(yNPArray)): #TODO CLEAN IT!
                    if xMin3 <= xNPArray[i] and xNPArray[i] <= xMax3 and yNPArray[i]>0.0 : 
                        tempYCenterList.append(yNPArray[i])
                yMinList3.append(min(tempYCenterList))
                yMaxList3.append(max(tempYCenterList))


    dx=x[1]-x[0]

    yMin = min(yMinList)/2
    yMax = max(yMaxList)*2
    xMin = min(xMinList)- 10 * dx
    xMax = max(xMaxList)+ 10 * dx

    #yMin2, Ymax2, xMin2, xMax2 = yMin, yMax, xMin, xMax
    
    limitsLabel = 'Signal limits'
    
    
    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2 ] = getHardcodedLimits(histName[0], wavelength, xMin, xMax)        
    ax1.set_xlabel(xLabel)
    #fig.text(0.5, 0.04, xLabel, ha='center')

    ax1.set_ylabel('Counts')
    
    
    if plotFit!=0:
        ax1.semilogy(x, gauss_function(x, *popt), '-', label='fitted curve', zorder=10)
         
    if plotLimits:
        if plotFit!=0 and gaussian==1: #else the limits are already set earlier
            limitLeft = meanFit - 3 * sigmaFit
            limitRight = meanFit + 3 * sigmaFit
        
        ax1.plot([limitLeft, limitLeft], [yMin, yMax], colors['magenta'])
        ax1.plot([limitRight, limitRight], [yMin, yMax], colors['magenta'], label=limitsLabel)
 

    ax1.set_ylim([yMin, yMax]) 
    ax1.set_xlim([xMin, xMax])
    ax1.set_yscale('log')
    ax1.legend(frameon=False)
    ax1.grid()

    if plotCenter:
        yMin2 = min(yMinList2)/2
        yMax2 = max(yMaxList2)*2 
            
        ax2.set_xlabel(xLabel)
        
        ax2.set_ylim([yMin2, yMax2]) 
        ax2.set_xlim([xMin2, xMax2])
        ax2.set_yscale('log')
        ax2.grid()

        ######################################
        yMin3 = min(yMinList3)/2
        yMax3 = max(yMaxList3)*2 
            
        ax3.set_xlabel(xLabel)
       
        ax3.set_ylim([yMin3, yMax3]) 
        ax3.set_xlim([xMin3, xMax3])
        ax3.set_yscale('log')
        ax3.grid()
    
    fig.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)
    plt.subplots_adjust(wspace=0.15)
    
    if createPDF:
        fig.savefig('gaussianFit_'+histName[0]+'.pdf')
    else:
        fig.savefig('gaussianFit_'+histName[0]+'.png')
    fig.show()
    
##########################################################################################################################################################
def fitGaussianTWO(shistFileName, histName, createPDF=0, normalise=1, plotFit=1, plotFitAll=0, plotLimits=1, plotCenter=0,):
    """Thisis just some hardcoded piece of ...code """

    setPlottingParameters(0.8)
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    if 1 < len(shistFileName):
        #plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']])))
        plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for material plots
        plt.rc('axes', prop_cycle=(cycler('color', [colors['blue'],  colors['darkred']]))) #for PE plot 0.6Ang
    else:
        plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['blue'], colors['black']])))
    plt.enable_tex_fonts()
    
    if plotCenter:
        fig, (ax1, ax2, ax3) = plt.subplots(1,3, gridspec_kw = {'width_ratios':[2,1,1]}, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
        
    else: 
        fig, ax1 = plt.subplots(1,1 , figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')

    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2] = getHardcodedLimits(histName[0], 11.0, 0, 0) #just to get xMin2, xMax2 !!
    xMin3, xMax3 = -0.01, 0.01

    xMinList = list()
    xMaxList = list()
    yMinList = list()
    yMaxList = list()
    yMinList2 = list()
    yMaxList2 = list()
    yMinList3 = list()
    yMaxList3 = list()
    legendList = list()

    tempVar = 0

    for fileIndex in range(len(shistFileName)):
        for histIndex in range(len(histName)):
            hc = sh.HistCollection(shistFileName[fileIndex])
            my_histogram = hc.hist(histName[histIndex])
    
            [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)

            x = np.array(histAxis)
            yErr = np.array(histError)
            normFactor = 1
            if normalise:
                normFactor = sum(histData)*(x[2]-x[1])
            y = [ tmp/normFactor for tmp in histData ]

            
            [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])

            
            if tempVar == 0 :
                labelText = 'Detection'
                tempVar = 1
            elif tempVar == 1 :
                labelText = 'No smearing'
                tempVar = 2
            else:
                labelText = 'Conversion'
                
            ax1.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText )

            yNPArray = np.array(y)
            xNPArray = np.array(x)
            
            yMinList.append(min(yNPArray[np.nonzero(yNPArray)])) 
            yMaxList.append(max(y))
            
            xMinList.append(min(xNPArray[np.nonzero(yNPArray)]))
            xMaxList.append(max(xNPArray[np.nonzero(yNPArray)]))
            
            if plotCenter:
                ax2.errorbar(x, y, yerr=yErr, fmt='o-', label=labelText , markersize=10)
                tempYCenterList = list()
                for i in range(len(yNPArray)): #TODO CLEAN IT!
                    if xMin2 <= xNPArray[i] and xNPArray[i] <= xMax2 and yNPArray[i]>0.0 : 
                        tempYCenterList.append(yNPArray[i])
                yMinList2.append(min(tempYCenterList))
                yMaxList2.append(max(tempYCenterList))


        
            if plotFit==1 or plotFitAll:
                #estimate mean and standard deviation
                mean = 0 # np.dot(x,y)
                sigma = 0.1 #this works for almost all 
                #sigma = 0.001# USE this if the fitting fails
                #do the fit!
                popt, pcov = curve_fit(gauss_function, x, y, p0 = [1, mean, sigma])
                ### plot fitted curve ###
                xFit = x
                #plt.semilogy(x, gauss_function(x, *popt), 'b-', label='fitted curve')
                meanFit = popt[1]
                sigmaFit = popt[2]

                print(shistFileName[fileIndex], histName[histIndex],'Mean: ', popt[1] ,  'Sigma ' ,  popt[2])
                plotFit=42 #Plot fitted curve only for the first

    
    for fileIndex in range(len(shistFileName)):
        for histogramName in ['neutron_dx_conv_center']:
            hc = sh.HistCollection(shistFileName[fileIndex])
            my_histogram = hc.hist(histogramName)
    
            [histData, histError, histAxis, titleText, histLabel] = disassembleHist1d(my_histogram)

            x = np.array(histAxis)
            yErr = np.array(histError)
            normFactor = 1
            if normalise:
                normFactor = sum(histData)*(x[2]-x[1])
            y = [ tmp/normFactor for tmp in histData ]

            
            [wavelength, panelNumber, material, ommit] = decomposeFileName(shistFileName[fileIndex])


            yNPArray = np.array(y)
            xNPArray = np.array(x)
            

            if plotCenter:
                ax3.errorbar(x, y, yerr=yErr, color=colors['black'], fmt='o-', label=labelText , markersize=10)
                tempYCenterList3 = list()
                for i in range(len(yNPArray)): #TODO CLEAN IT!
                    if xMin3 <= xNPArray[i] and xNPArray[i] <= xMax3 and yNPArray[i]>0.0 : 
                        tempYCenterList3.append(yNPArray[i])
                yMinList3.append(min(tempYCenterList3))
                yMaxList3.append(max(tempYCenterList3))

                
    dx=x[1]-x[0]

    yMin = min(yMinList)/2
    yMax = max(yMaxList)*2
    xMin = min(xMinList)- 10 * dx
    xMax = max(xMaxList)+ 10 * dx

    #yMin2, Ymax2, xMin2, xMax2 = yMin, yMax, xMin, xMax
    
    limitsLabel = 'Signal limits'
    
    
    [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2 ] = getHardcodedLimits(histName[0], wavelength, xMin, xMax)        
    ax1.set_xlabel(xLabel)
    #fig.text(0.5, 0.04, xLabel, ha='center')

    ax1.set_ylabel('Counts')
    
    
    if plotFit!=0:
        ax1.semilogy(x, gauss_function(x, *popt), '-', label='fitted curve', zorder=10)
         
    if plotLimits:
        if plotFit!=0 and gaussian==1: #else the limits are already set earlier
            limitLeft = meanFit - 3 * sigmaFit
            limitRight = meanFit + 3 * sigmaFit
        
        ax1.plot([limitLeft, limitLeft], [yMin, yMax], colors['magenta'])
        ax1.plot([limitRight, limitRight], [yMin, yMax], colors['magenta'], label=limitsLabel)
 

    ax1.set_ylim([yMin, yMax]) 
    ax1.set_xlim([xMin, xMax])
    ax1.set_yscale('log')
    ax1.legend(frameon=False, loc='upper left')
    ax1.grid()

    if plotCenter:
        yMin2 = min(yMinList2)/2
        yMax2 = max(yMaxList2)*2 
        if plotFit!=0:
            ax2.semilogy(x, gauss_function(x, *popt), 'b-', label='fitted curve')
        if plotLimits:
            ax2.plot([limitLeft, limitLeft], [yMin2, yMax2], colors['magenta'])
            ax2.plot([limitRight, limitRight], [yMin2, yMax2], colors['magenta'], label=limitsLabel)
            
        ax2.set_xlabel(xLabel)
        #ax2.yaxis.tick_right() #is it better?
        
        ax2.set_ylim([yMin2, yMax2]) 
        ax2.set_xlim([xMin2, xMax2])
        ax2.set_yscale('log')
        ax2.grid()

        yMin3 = min(yMinList3)/2
        yMax3 = max(yMaxList3)*2 
        if plotFit!=0:
            ax3.semilogy(x, gauss_function(x, *popt), 'b-', label='fitted curve')
        if plotLimits:
            ax3.plot([limitLeft, limitLeft], [yMin3, yMax3], colors['magenta'])
            ax3.plot([limitRight, limitRight], [yMin3, yMax3], colors['magenta'], label=limitsLabel)
            
        ax3.set_xlabel(xLabel)
        #ax3.yaxis.tick_right() #is it better?
        
        ax3.set_ylim([yMin3, yMax3]) 
        ax3.set_xlim([xMin3, xMax3])
        ax3.set_yscale('log')
        ax3.grid()
    
    fig.tight_layout(pad=0.1, w_pad=0.2, h_pad=0.1)
    plt.subplots_adjust(wspace=0.15)
    
    if createPDF:
        fig.savefig('gaussianFit_'+histName[0]+'.pdf')
    else:
        fig.savefig('gaussianFit_'+histName[0]+'.png')
    fig.show()

##########################################################################################################################################################
def getHardcodedLimits(histName, wavelength, xMin2, xMax2):

    limitLeft, limitRight = 0, 0
    gaussian=1
    
    if 'dtof' in histName: #lambda dependent
        gaussian=0
        xLabel = '{$\delta$}TOF [ms]'
        lambdaLimit ={ #calculated with s=5m
            "0.6": 0.000565,
            "1.8": 0.001695,
            "3.0": 0.002825,
            "5.0": 0.004708,
            "11.0":0.01063
        }
        limitLeft, limitRight = -lambdaLimit[str(wavelength)], lambdaLimit[str(wavelength)]
        xMin2, xMax2 = -0.02, 0.02
    if 'dlambda' in histName: #lambda dependent
        gaussian=0
        xLabel = '{$\delta\lambda$} [{\AA}]'

        lambdaLimit ={ #calculated with s=5m
            "0.6": 0.000447,
            "1.8": 0.001341,
            "3.0": 0.002235,
            "5.0": 0.003725,
            "11.0":0.008195
        }
        limitLeft, limitRight = -lambdaLimit[str(wavelength)], lambdaLimit[str(wavelength)]
        xMin2, xMax2 = -0.013, 0.013
        #xMin2, xMax2 = -0.005, 0.005
    if 'dekin' in histName: #lambda dependent
        gaussian=0
        xLabel = '{$\delta$E} [meV]'

        lambdaLimit ={ #NOPE NOT REAL DATA
            "0.6": 0.000447,
            "1.8": 0.001341,
            "3.0": 0.002235,
            "5.0": 0.003725,
            "11.0":0.008195
        }
        limitLeft, limitRight = -lambdaLimit[str(wavelength)], lambdaLimit[str(wavelength)]
        xMin2, xMax2 = -0.5, 1.0
    elif 'dth' in histName:
        limitLeft, limitRight = 3* -0.0345141928622, 3* 0.0345141928622
        xLabel = '{$\delta\Theta$} [degree]'
        xMin2, xMax2 = -0.3, 0.3
    elif 'dphi' in histName:
        xLabel = '{$\delta\Phi$} [degree]'
        limitLeft, limitRight = 3* -0.542488871037, 3* 0.542488871037
        xMin2, xMax2 = -10, 10
    elif 'dx' in histName:
        xLabel = '{$\delta$}X [cm]'
        limitLeft, limitRight = 3* -0.291508708298, 3* 0.291508708298
        xMin2, xMax2 =-2, 2
            
    elif 'dQ' in histName:
        xLabel = '{$\delta$}Q [1/{\AA}]'
        # limitLeft, limitRight = 3* -0.00210002220634,3* 0.00210002220634
        
        xMin2, xMax2 = -0.02, 0.02
    elif 'dy' in histName:
        gaussian=0
        xLabel = '{$\delta$}Y [cm]'
        limitLeft, limitRight = -0.3725, 0.3725
        xMin2, xMax2 = -1, 1

    return [limitLeft, limitRight, gaussian, xLabel, xMin2, xMax2]
    
    
    ###########################################################################################################################################################
def materialAbsorptionOverlay(fileNames, createPDF=0):
    valuesAl = list()
    valuesCu = list()
    valuesAlCu = list()
    valuesB = list()
    errorAl = list()
    errorCu = list()
    errorAlCu = list()
    errorB = list()
    wavelengthAl = list()
    wavelengthCu = list()
    wavelengthAlCu = list()
    wavelengthB = list()

    def calcError(x):
        N=2e7
        return m.sqrt(x*100/N) 
 
    for i in range(len(fileNames)):
        [wavelength, panelNumber, material, ommit] = decomposeFileName(fileNames[i])
        
        hc = sh.HistCollection(fileNames[i])
        tubeWall_hit = float(hc.hist('h_st_tubeWall_hit').getIntegral())
        strawWall_hit = float(hc.hist('h_st_strawWall_hit').getIntegral())
        sumAbsorption = tubeWall_hit + strawWall_hit
        numberOfNeutrons = float(hc.hist('neutron_true_x').getIntegral())
        #numberOfHits = float(hc.hist('neutron_x_hit').getIntegral())
        relativeAbsorption = sumAbsorption/numberOfNeutrons

        print(fileNames[i]+' Al: '+ str(tubeWall_hit))
        print(fileNames[i]+' Cu: '+ str(strawWall_hit))
        

        if material == 'Al':
            valuesAl.append(relativeAbsorption)
            errorAl.append(calcError(relativeAbsorption))
            wavelengthAl.append(wavelength)
        if material == 'Cu':
            valuesCu.append(relativeAbsorption)
            errorCu.append(calcError(relativeAbsorption))
            wavelengthCu.append(wavelength)
        if material == 'AlCu':
            valuesAlCu.append(relativeAbsorption)
            errorAlCu.append(calcError(relativeAbsorption))
            wavelengthAlCu.append(wavelength)
        if material == 'B':
            numberOfConversions = float(hc.hist('neutron_x_conv').getIntegral())
            numberOfHits = float(hc.hist('neutron_x_hit').getIntegral())
            boronAbsorption = numberOfConversions - numberOfHits
            relativeAbsorption = boronAbsorption / numberOfNeutrons
            #print(relativeAbsorption)
            valuesB.append(relativeAbsorption)
            errorB.append(calcError(relativeAbsorption))
            wavelengthB.append(wavelength)
        
        numberOfConversions = float(hc.hist('neutron_x_conv').getIntegral())
        numberOfHits = float(hc.hist('neutron_x_hit').getIntegral())
        boronAbsorption = numberOfConversions - numberOfHits
        print('number of conv: ' + str(numberOfConversions))
        print('number of hits: ' + str(numberOfHits))
        print(fileNames[i]+' B: '+ str(boronAbsorption))
        print('DetToConv ration: '+ str(numberOfHits/numberOfConversions))
           
        print(' ')
        

    setPlottingParameters(multiplicationFactor = 0.7)
    #setPlottingParameters()
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']])))

    #plt.errorbar(xAxisValues, values0p6, yerr=errors0p6, fmt='o--')
    
    #plt.plot(wavelengthB, valuesB, '--o', linewidth=2, markersize=6)
    #plt.plot(wavelengthAlCu, valuesAlCu, '--o', linewidth=2, markersize=6)
    #plt.plot(wavelengthCu, valuesCu, '--o', linewidth=2, markersize=6)
    #plt.plot(wavelengthAl, valuesAl, '--o', linewidth=2, markersize=6)

    plt.errorbar(wavelengthB, valuesB, yerr=errorB, fmt='o--', linewidth=2, markersize=6)
    plt.errorbar(wavelengthAlCu, valuesAlCu, yerr=errorAlCu, fmt='o--', linewidth=2, markersize=6)
    plt.errorbar(wavelengthCu, valuesCu, yerr=errorCu, fmt='o--', linewidth=2, markersize=6)
    plt.errorbar(wavelengthAl, valuesAl, yerr=errorAl, fmt='o--', linewidth=2, markersize=6)
    
    #plt.legend(['Al+Cu','Cu', 'Al'], frameon=False)
    #plt.axis([0, 12, 0, 1.1*max([max(valuesAl),max(valuesCu),max(valuesAlCu)]) ])
    plt.legend([ 'B', 'Al+Cu','Cu', 'Al'], frameon=False)
    plt.axis([0, 12, 0, 1.1*max([max(valuesAl),max(valuesCu),max(valuesAlCu),max(valuesB)]) ])
    plt.xlabel('Wavelength [{\AA}]')
    plt.ylabel('Relative absorption')
    plt.grid()
    plt.tight_layout(pad=0.28, w_pad=0.2, h_pad=0.1)
    
    #plt.show()

    if createPDF:
        plt.savefig('materialAbsorption.pdf')
    else:
        plt.savefig('materialAbsorption.png')

    plt.show()


##########################################################################################################################################################
def plotCompareHardCodedPanelEfficiency(createPDF=1):

    xAxisQantity = 'panelNumber'
    
    values0p6 = [6.50, 12.61, 18.11, 23.00, 27.27]
    values1p8 = [15.24, 27.74, 36.98, 44.17, 49.56]
    values3 = [22.95, 39.01, 49.06, 55.41, 59.59]
    values5 = [31.36, 49.69, 58.41, 62.75, 65.04]
    values11 = [45.57, 62.07, 65.34, 66.16, 66.37]

    def calcError(x):
        N=2e7
        return m.sqrt(x*100/N) 
    
    errors0p6 = [calcError(6.50), calcError(12.61), calcError(18.11), calcError(23.00), calcError(27.27)]
    errors1p8 = [calcError(15.24), calcError(27.74), calcError(36.98), calcError(44.17), calcError(49.56)]
    errors3 = [calcError(22.95), calcError(39.01), calcError(49.06), calcError(55.41), calcError(59.59)]
    errors5 = [calcError(31.36), calcError(49.69), calcError(58.41), calcError(62.75), calcError(65.04)]
    errors11 = [calcError(45.57), calcError(62.07), calcError(65.34), calcError(66.16), calcError(66.37)]
    
    
    xAxisValues = [1,2,3,4,5]
    
    legendList = ['$\lambda$=0.6 {\AA}','$\lambda$=1.8 {\AA}', '$\lambda$=3 {\AA}', '$\lambda$=5 {\AA}', '$\lambda$=11 {\AA}']
    minValueList = list()
    maxValueList = list()

    fig = plt.figure()
    #fig.set_size_inches(12,10.2)
    #setPlottingParameters(1)
    #plt.rc('legend',**{'fontsize':30})
    #plt.rc('legend',**{'labelspacing':0.3})
    fig.set_size_inches(8, 6)
    setPlottingParameters(0.8)
    plt.tight_layout()
    
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['cyan'], colors['navy']]))) #for 5 wavelength
    
    mSize=10
    lWidth=4

    plt.errorbar(xAxisValues, values0p6, yerr=errors0p6, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values1p8, yerr=errors1p8, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values3, yerr=errors3, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values5, yerr=errors5, fmt='o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(xAxisValues, values11, yerr=errors11, fmt='o--', linewidth=lWidth, markersize=mSize)
       
    plt.legend(legendList, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.03))
    plt.ylabel('Detection efficiency [\%]')
    plt.grid()

    yAxisMin = 0
    yAxisMax = 100
    xAxisMin = 0.5
    xAxisMax = 5.5
        
    plt.xlabel('Panel number')
    plt.axis([xAxisMin, xAxisMax, yAxisMin, yAxisMax])

    plt.plot([xAxisMin, xAxisMax], [70, 70], colors['magenta'], linewidth=4) #, label=limitsLabel)
    
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        plt.savefig('efficiencyPanel.pdf')
    else:
        plt.savefig('efficiencyPanel.png')
    plt.show()


def pulseHeightSpectra(fileName, createPDF=0):

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color',  [colors['orangered'], colors['darkred'], colors['cyan'], colors['navy']]))) #,  colors['darkgreen'],
    #plt.rc('axes', prop_cycle=(cycler('color',  [colors['orangered'], colors['darkred'], colors['navy'], colors['darkgreen']]))) #, colors['cyan'],

    hc = sh.HistCollection(fileName)
    setPlottingParameters(0.7)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.tight_layout()

    mSize=6
    lWidth=3

    histNames = [ 'energyDeposition_lithium', 'energyDeposition_alpha','energyDeposition_alphalithium','energyDeposition']
    legendNames = ['$^7$Li',r'$\alpha$',r'$^7$Li + $\alpha$',r'$^7$Li + $\alpha$ + $\gamma$',]
    #histNames = [ 'energyDeposition_lithium', 'energyDeposition_alpha','energyDeposition_alphalithium']
    #legendNames = ['$^7$Li',r'$\alpha$',r'$^7$Li + $\alpha$']
  
    for histname in histNames:

        hist1D = hc.hist(histname)
        [histData, histError, histAxis, histTitleText, histLabel] = BCSplt.disassembleHist1d(hist1D)
        ax.errorbar(histAxis, histData, yerr=histError, fmt='o--', linewidth=lWidth, markersize=mSize)


    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend(legendNames, frameon=False)
    ax.set_ylabel('Counts')
    ax.grid()

    #ax.set_xlabel('Energy [keV]')
    ax.set_xlabel('Energy [MeV]')
    ax.set_xlim(-30, 2000)
    ax.set_xticks([0, 400, 800, 1200, 1600, 2000 ])
    #ax.set_xticklabels(["0","400","800","1200","1600","2000"])
    ax.set_xticklabels(["0","0.4","0.8","1.2","1.6","2.0"])
    #ax.set_xticks([0, 400, 800, 1200, 1600, 2000, 2400, 2800 ])
    #ax.set_xticklabels(["0","400","800","1200","1600","2000","2400","2800"])
    ylimMin = -5000
    ylimMax = 125000
    ax.set_ylim(ylimMin, ylimMax)


    axTop = ax.twiny()
    axTop.xaxis.set_ticks_position('top')
    axTop.set_xlim(-30, 2000)
    axTop.set_xticks([840, 1010, 1470, 1780 ])
    #axTop.set_xticklabels(["840","1010","1470","1780"])
    axTop.set_xticklabels(["0.84","1.01","1.47","1.78"])
    axTop.set_ylim(ylimMin, ylimMax)
    
    axTop.plot([840, 840], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1010, 1010], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1470, 1470], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1780, 1780], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    
    fig.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    setPlottingParameters(0.8)
    
    plt.savefig('pulseHeightSpectra.pdf')
    #plt.show()
    
def pulseHeightSpectraLog(fileName, createPDF=0):

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    #plt.rc('axes', prop_cycle=(cycler('color',  [colors['orangered'], colors['darkred'], colors['cyan'], colors['navy']]))) #,  colors['darkgreen'],
    plt.rc('axes', prop_cycle=(cycler('color',  [colors['magenta'], colors['orangered'], colors['darkred'], colors['cyan'], colors['navy']]))) #,  
    
    hc = sh.HistCollection(fileName)
    setPlottingParameters(0.7)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.tight_layout()

    mSize=6
    lWidth=3
    
    histNames = [ 'energyDeposition_gammas', 'energyDeposition_lithium', 'energyDeposition_alpha','energyDeposition_alphalithium','energyDeposition']
    legendNames = ['$\gamma$','$^7$Li',r'$\alpha$',r'$^7$Li + $\alpha$',r'$^7$Li + $\alpha$ + $\gamma$',]
  
    for histname in histNames:

        hist1D = hc.hist(histname)
        [histData, histError, histAxis, histTitleText, histLabel] = BCSplt.disassembleHist1d(hist1D)
        ax.errorbar(histAxis, histData, yerr=histError, fmt='o--', linewidth=lWidth, markersize=mSize)


    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend(legendNames, frameon=False)
    ax.set_ylabel('Counts')
    ax.grid()

    #ax.set_xlabel('Energy [keV]')
    ax.set_xlabel('Energy [keV]')
    ax.set_xlim(-30, 2800)
    ax.set_xticks([0, 400, 800, 1200, 1600, 2000, 2400, 2800 ])
    #ax.set_xticklabels(["0","400","800","1200","1600","2000","2400","2800"])
    ax.set_xticklabels(["0","0.4","0.8","1.2","1.6","2.0","2.4","2.8"])

    ylimMax = 1e7
    ylimMin = 0.1
    ax.set_ylim(ylimMin, ylimMax)
    ax.set_yscale('log') 

    axTop = ax.twiny()
    axTop.xaxis.set_ticks_position('top')
    axTop.set_xlim(-30, 2800)
    axTop.set_xticks([840, 1010, 1470, 1780 ])
    #axTop.set_xticklabels(["840","1010","1470","1780"])
    axTop.set_xticklabels(["0.84","1.01","1.47","1.78"])
    axTop.set_yscale('log') 
    axTop.set_ylim(ylimMin, ylimMax)
    
    axTop.plot([840, 840], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1010, 1010], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1470, 1470], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    axTop.plot([1780, 1780], [ylimMin, ylimMax], 'k--', alpha=0.5, linewidth=1 )
    
    fig.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    setPlottingParameters(0.8)
    
    plt.savefig('pulseHeightSpectraLog.pdf')
    #plt.show()
    

    ##########################################################################################################################################################
def plotHardCodedEfficiencyValidation(createPDF=1):


    measurementDataMatrix = numpy.array([[1.5586760280842529, 44.215448964122466],
                                         [1.6549648946840518, 42.38443793758593],
                                         [2.3289869608826477, 55.73966173329148],
                                         [2.6178535606820463, 56.17901687614643],
                                         [3.0992978936810434, 58.190810162220565],
                                         [3.709127382146439, 59.6805574848865],
                                         [4.968906720160481, 69.98873392259657],
                                         [5.739217652958876, 70.69506228609492],
                                         [6.61384152457372, 70.35563615601438],
                                         [9.037111334002006, 73.26101532514666],
                                         [10.369107321965897, 72.22864886918116],
                                         [11.580742226680039, 72.59082591285308]])


    mcnpDataMatrix = numpy.array([[1.0050150451354058, 30.163609694951816],
                                  [1.5105315947843527, 39.76563277619117],
                                  [1.7993981945837505, 44.39255607607993],
                                  [2.0080240722166494, 47.18654327212827],
                                  [2.5055165496489464, 52.775392699363074],
                                  [3.0110330992978933, 57.230196587581716],
                                  [4.006018054162487, 63.52239925884523],
                                  [5.00100300902708, 67.3718538385057],
                                  [6.004012036108325, 69.91278087260692],
                                  [7.007021063189569, 71.58129787399275],
                                  [8.00200601805416, 72.72628135223556],
                                  [9.005015045135405, 73.43514731763447],
                                  [10.00802407221665, 73.96953127649033],
                                  [10.994984954864593, 74.1547762152324],
                                  [12.006018054162485, 74.42752466778745]])


    davideDataMatrix = numpy.array([[0.10892018779342733, 0.048858478797055516],
                                    [0.34553990610328644, 0.11452370749173879],
                                    [0.5788732394366198, 0.1655852355512386],
                                    [0.815492957746479, 0.22313687398223025],
                                    [1.0521126760563382, 0.2668954089649459],
                                    [1.2887323943661975, 0.3155220981058767],
                                    [1.5253521126760563, 0.35116704282490085],
                                    [1.7586854460093897, 0.37951051814606374],
                                    [1.9953051643192485, 0.42002361702330293],
                                    [2.225352112676056, 0.4524246493157731],
                                    [2.4619718309859153, 0.44587892466360035],
                                    [2.6953051643192487, 0.4815246312220859],
                                    [2.928638497652582, 0.48958413088401953],
                                    [3.1718309859154927, 0.5219821158186442],
                                    [3.4051643192488257, 0.5032667676103953],
                                    [3.638497652582159, 0.5194398575360208],
                                    [3.875117370892019, 0.48125113085545057],
                                    [4.111737089201878, 0.5436709234446572]])

    
    def calcError(x):
        N=2e7
        return m.sqrt(x*100/N) 
    
    ownSimulationDataMatrix = numpy.array([[0.6, 27.27],[1.8, 49.56],[3, 59.59],[5, 65.04],[11, 66.37]])
    ownSimulationError = numpy.array([calcError(27.27), calcError(49.56), calcError(59.59), calcError(65.04), calcError(66.37)])


    
    #legendList = ['Measurement','MCNP simulation', 'Geant4 simulation']
    legendList = ['Measurement 1','MCNP simulation', 'Measurement 2','Geant4 simulation']
    minValueList = list()
    maxValueList = list()

    fig = plt.figure()
    #fig.set_size_inches(12,10.2)
    #setPlottingParameters(1)
    #plt.rc('legend',**{'fontsize':30})
    #plt.rc('legend',**{'labelspacing':0.3})
    fig.set_size_inches(8, 6)
    setPlottingParameters(0.8)
    plt.tight_layout()
    
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    plt.rc('axes', prop_cycle=(cycler('color', [colors['darkred'], colors['orangered'],  colors['darkgreen'], colors['navy'], colors['cyan']]))) #for 5 wavelength
    
    mSize=9
    lWidth=3


    plt.plot(measurementDataMatrix[:,0], measurementDataMatrix[:,1], 'o--', linewidth=lWidth, markersize=mSize)
    plt.plot(mcnpDataMatrix[:,0], mcnpDataMatrix[:,1], 'o--', linewidth=lWidth, markersize=mSize)
    plt.plot(davideDataMatrix[:,0], davideDataMatrix[:,1]*100, 'o--', linewidth=lWidth, markersize=mSize)
    plt.errorbar(ownSimulationDataMatrix[:,0], ownSimulationDataMatrix[:,1], yerr=ownSimulationError, fmt='o--', linewidth=lWidth, markersize=mSize)
    
                                 
    #plt.errorbar(xAxisValues, values1p8, yerr=errors1p8, fmt='o--', linewidth=lWidth, markersize=mSize)
    
       
    plt.legend(legendList, frameon=False, loc='lower right') #bbox_to_anchor=(0, 1.03))
    plt.ylabel('Detection efficiency [\%]')
    plt.grid()

    yAxisMin = 0
    yAxisMax = 100
    xAxisMin = 0
    xAxisMax = 12.5


    plt.yticks(np.arange(0, 105, step=10))
    plt.xticks(np.arange(0, 12.5, step=2)) 
    
    plt.xlabel('Neutron wavelength [{\AA}]')
    plt.axis([xAxisMin, xAxisMax, yAxisMin, yAxisMax])

    #plt.plot([xAxisMin, xAxisMax], [70, 70], colors['magenta'], linewidth=4) #, label=limitsLabel)
    
    plt.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.1)

    if createPDF:
        plt.savefig('efficiencyValidation.pdf')
    else:
        plt.savefig('efficiencyValidation.png')
    plt.show()

 

