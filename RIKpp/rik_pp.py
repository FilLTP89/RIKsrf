# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Set of subroutines to post-process RIK output
"""

#=======================================================================
# Required modules
#=======================================================================
# System/OS
import os
import sys
sys.path.append("/home/filippo/Data/Filippo/RIKsrf/RIKpp")
from sys import exit
# SciPy/NumPy
import pylab as p
import numpy as np
import scipy.interpolate
from scipy import signal
from scipy.interpolate import interp1d
from INTFORTRAN import trapezoid
import math
import shutil
import h5py

# Plotting
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib import ticker
from matplotlib.ticker import LogFormatter
from matplotlib.colors import LogNorm,LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon,Circle
from matplotlib.collections import PatchCollection
from matplotlib.image import NonUniformImage
from collections import OrderedDict

#=======================================================================
# General informations
#=======================================================================
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2018, MSSMat UMR CNRS 8579 - CentraleSupélec"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"


def moment_computation(M0, time, f, ts, gamma):

    T = 1.0 / f
    moment = np.zeros(len(time))
    for i in np.arange(len(time)):
        t = time[i]
        if (t >= ts):
            s = ((t - ts) / T)**gamma
            moment[i] = M0 * (1.0 - (1.0 + s) * math.e**(-s))

    return moment

def yeni_data_seti(data, dim, chunk, attrib, output):
    ''' Attribuer un set de data de dimension dim
    comme attribut de attrib dans un fichier output de hdf5 '''
    dataset = output.create_dataset(attrib, dim, chunks=chunk)
    dataset = data

def vecteurs(aD, aS, aR, output):
    ''' Calcul du vecteur normal Vnormal sur le plan de faille, et du vecteur unitaire
    de glissement Vslip dans le repere de reference et les attribuer dans le fichier 
    output de hdf5 '''

    # En supposant que (puisque l'on multiplie avec MatMesh)
    # x ----> EST
    # y ----> NORD
    # z ----> UP

    Vnormal = np.array([+np.sin(aD) * np.cos(aS),
                        -np.sin(aD) * np.sin(aS),
                        +np.cos(aD)])

    Vslip   = np.array([-np.sin(aR) * np.cos(aD) * np.cos(aS) + 
                        +np.cos(aR) * np.sin(aS),
                        +np.sin(aR) * np.cos(aD) * np.sin(aS) +
                        +np.cos(aR) * np.cos(aS),
                        +np.sin(aR) * np.sin(aD)])

    output.attrs['Vnormal'] = Vnormal
    output.attrs['Vslip'] = Vslip
    print "Vector normal to the fault : ", Vnormal
    print "Vector of the slip         : ", Vslip

    ###
    M = np.zeros((3, 3))
    for i in np.arange(3):
        for j in np.arange(3):
            M[i, j] = Vnormal[i] * Vslip[j] + Vnormal[j] * Vslip[i]
    print 'Moment matrix of SEM3D code: '
    print M[0, :]
    print M[1, :]
    print M[2, :]

def read_moment_rate_RIK(dosya, dim):
    ''' Lire la vitesse de moment par un fichier dosya et 
    l'assigner a l'array de dimension dim '''

    Mrate = np.zeros((dim[0] * dim[1], dim[2]))

    print 'Reading moment rate of all points...'
    f = open(dosya, 'r')
    # Pour tous les points
    for i in np.arange(dim[0] * dim[1]):
        # Histoire temporelle
        for t in np.arange(dim[2]):
            string = f.readline()
            dagit = string.rsplit()
            Mrate[i, t] = float(dagit[1])

            f.readline()
            f.readline()

    Mratenew = np.reshape(Mrate, dim, order='F')
    return Mratenew

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def readhypo(filename):
    x = np.genfromtxt(filename, usecols=0)
    y = np.genfromtxt(filename, usecols=1)
    return x, y

def read_subsources(filename):
    f = open(filename, 'r')
    print 'Reading the file', filename
    print '...'
    xgrid = np.genfromtxt(filename, usecols=0)
    ygrid = np.genfromtxt(filename, usecols=1)
    ssources = np.genfromtxt(filename, usecols=2)

    return xgrid, ygrid, ssources

def readfileslip(NSR, filename):
    xgrid = np.zeros((NSR))
    ygrid = np.zeros((NSR))
    slip = np.zeros((NSR))

    f = open(filename, 'r')
    print 'Reading the file', filename
    print '...'

    for i in np.arange(NSR):
        string = f.readline()
        degerler = [float(j) for j in string.split()]
        xgrid[i] = degerler[0]
        ygrid[i] = degerler[1]
        slip[i] = degerler[2]
    return xgrid, ygrid, slip

def plotslip(nrows, ncols, LF, WF, slip, kaydet, figname, nucx, nucy):
    print 'Plotting max-slip distribution'
    print
    # Making a colormap
    c = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
          '#fcab6c', '#f50000']
    #cc1 = np.logspace(np.log10(0.25),np.log10(0.95),6)
    cc1 = np.linspace(0, 1, 6)
    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0],
                          c(cc[1]), c(cc[2]), cc1[1],
                          c(cc[2]), c(cc[3]), cc1[2],
                          c(cc[3]), c(cc[4]), cc1[3],
                          c(cc[4]), c(cc[5]), cc1[4],
                          c(cc[5]), c(cc[6]), cc1[5],
                          c(cc[6])])
    # Figure parameters
    sns.set_style('whitegrid')
    fig = p.figure(figsize=(12,6))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Along strike [km]', fontsize=17)
    ax.set_ylabel('Along up-dip [km]', fontsize=17)
    print LF,WF
    ax.set_xlim([0.0,LF])
    ax.set_ylim([0.0,WF])
    vmin = slip.min()
    vmax = slip.max()
    print 'min and max of slip:'
    print vmin, vmax
    print

    grid = slip.reshape((nrows, ncols))
    im = plt.imshow(grid, extent=(0.0,LF,0.0,WF),
                    interpolation='bilinear', cmap=cmap, origin='lower')
    formatter = LogFormatter(10, labelOnlyBase=False)
    cb = fig.colorbar(im, shrink=0.5, aspect=10, pad=0.01,
                      ticks=np.arange(vmin, vmax, 1), format=formatter)
    cb.set_label('Slip [m]', labelpad=20, y=0.5, rotation=90, fontsize=17)
    p.setp(cb.ax.yaxis.get_ticklabels(), fontsize=16)
    ax.plot(nucx, nucy, marker='*', color='red', markersize=20)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    topname = 'Max slip [m] = ' + '% .4f' % max(slip)
    plt.title(topname + '\n', fontsize=20)
    if kaydet:
        fig.savefig(figname,dpi=500,bbox_inches='tight')
    else:
        plt.show()

def plotsubsources(x, y, radius, LF, WF, kaydet, figname, nucx, nucy):

    print 'Plotting sub-sources'
    print
    # Making a colormap
    c = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
            '#fcab6c', '#f50000']
    #cc1 = np.logspace(np.log10(0.25),np.log10(0.95),6)
    cc1 = np.linspace(0, 1, 6)
    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0],
        c(cc[1]), c(cc[2]), cc1[1],
        c(cc[2]), c(cc[3]), cc1[2],
        c(cc[3]), c(cc[4]), cc1[3],
        c(cc[4]), c(cc[5]), cc1[4],
        c(cc[5]), c(cc[6]), cc1[5],
        c(cc[6])])
    # Figure parameters
    sns.set_style('whitegrid')
    fig = p.figure(figsize=(18, 10))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Along strike [km]', fontsize=17)
    ax.set_ylabel('Along up-dip [km]', fontsize=17)
    ax.set_xlim([0, 1.1 * LF])
    ax.set_ylim([-0.2, 1.1 * WF])
    i = 0
    for xx, yy in zip(x, y):
        circ = Circle((xx, yy), radius[i], fill=False, lw=1)
        ax.add_patch(circ)
        i = i + 1

    Nsub = len(radius)
    topname = 'Total number of subsources = ' + '%d' % Nsub
    plt.title(topname + '\n', fontsize=20)
    ax.plot(nucx, nucy, marker='*', color='red', markersize=20)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    if kaydet:
        fig.savefig(figname,dpi=500,bbox_inches='tight')
    else:
        plt.show()


def plotgrid(nrows, ncols, LF, WF, kaydet, figname, hypox, hypoy, x, y):

    print 'Creating the 2D grid...'
    # Making a colormap
    c = mcolors.ColorConverter().to_rgb
    cc = ['#ffffff', '#dfe6ff', '#a5b5da', '#516b8e', '#c5ce9b',
          '#fcab6c', '#f50000']
    #cc1 = np.logspace(np.log10(0.25),np.log10(0.95),6)
    cc1 = np.linspace(0, 1, 6)
    cmap = make_colormap([c(cc[0]), c(cc[1]), cc1[0],
                          c(cc[1]), c(cc[2]), cc1[1],
                          c(cc[2]), c(cc[3]), cc1[2],
                          c(cc[3]), c(cc[4]), cc1[3],
                          c(cc[4]), c(cc[5]), cc1[4],
                          c(cc[5]), c(cc[6]), cc1[5],
                          c(cc[6])])
    # Figure parameters
    sns.set_style('whitegrid')
    fig = p.figure(figsize=(18, 10))
    p.subplots_adjust(hspace=0.35)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Along strike [km]', fontsize=17)
    ax.set_ylabel('Along up-dip [km]', fontsize=17)
    ax.set_xlim([0.0,LF])
    ax.set_ylim([0.0,WF])

    plt.scatter(x, y, s=2, c='gray')
    ax.plot(hypox, hypoy, marker='*', color='red', markersize=25)

    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)

    topname = '2D grid for RIK model'
    plt.title(topname + '\n', fontsize=20)

    if kaydet:
        fig.savefig(figname,dpi=500,bbox_inches='tight')
    else:
        plt.show()