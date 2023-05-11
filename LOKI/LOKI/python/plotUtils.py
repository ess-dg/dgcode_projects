from __future__ import print_function

from builtins import map
from builtins import range
#matplotlib.use('agg')
##import PhysValTS.tns_style
from PyAna import *
#plt.enable_tex_fonts()

import numpy  as np

import operator
import matplotlib.pyplot as plt
from cycler import cycler

#from matplotlib import colors as mcolors #added lated, might cause problems
#from matplotlib.pyplot import figure

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

# def _col(cmap, i, n):
#     assert i < n
#     return cmap(i / (n - 1.0)) if n > 1 else cmap(1.0)
# cmap = _defcolmap

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

    #colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
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

def disassembleHist1d(hist1d):
    """This function returns the Data, Error, AxisX, Title and LabelX of a hist1D object in separate arrays/variables"""

    nBins = hist1d.getNBins()

    histData = [hist1d.getBinContent(i) for i in range(nBins)]
    histError = [hist1d.getBinError(i) for i in range(nBins)]
    histAxis = [hist1d.getBinCenter(i) for i in range(nBins)]
    titleText = hist1d.getTitle()
    label = hist1d.getXLabel()

    return (histData, histError, histAxis, titleText, label)

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
