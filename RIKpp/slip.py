# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Script to plot 3 subplots in a row
from /home/filippo/Data/Filippo/RIKsrf/RIKpp/rik_pp import *

if __name__ == '__main__':

    # Input data
    Lgrid   = 210
    Wgrid   = 170
    NSR = Wgrid*Lgrid

    filename = 'slipdistribution.dat'
    x, y, slip = readfile (NSR, filename)

    # Attention a l'ordre de W et L #
    plotslip (Wgrid,Lgrid,slip,False,'xxx.png',24.55,7.706,
            [0.0,29.4],[0.0,21.25])

    print '*** Done ***'
    #
