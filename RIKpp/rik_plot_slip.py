# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
Main file to plot slip patches obtained with RIKsrf2
'''
#=======================================================================
# Required modules
#=======================================================================
# System/OS
import os
import sys
sys.path.append("/home/filippo/Data/Filippo/RIKsrf/RIKpp")
from rik_pp import *

if __name__ == '__main__':

    # Input data
    Lgrid   = 294
    Wgrid   = 213
    NSR = Wgrid*Lgrid

    filename = 'slipdistribution.dat'
    x, y, slip = readfile (NSR, filename)

    # Attention a l'ordre de W et L #
    plotslip (Wgrid,Lgrid,slip,False,'xxx.png',24.55,7.7,
            [0.0,29.4],[0.0,21.25])

    print '*** Done ***'
