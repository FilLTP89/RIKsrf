# ********************************************************* #
#  Filename: 2d_plot_gallovics_output.py					#
#  Purpose : Script to plot 2D plots of Gallovic model      #
#  Use the instructions below 'Input data'                  #
#  Author  : Elif ORAL 										#
#  Email   : elif.oral@centralesupelec.fr 					#
# ********************************************************* #

# Librairies necessaires 
from elfolibsource import *



###############
### PROGRAM ###
###############


# Fichiers necessaires:
#######################
# slipdistribution.dat
# nucleationpoint.dat
# subsources.dat
#


# ********************************************************* #
# Input data
# ********************************************************* #

# Nombre de points de grille (toute la grille si different 
# que strong-motion generation area)
Lgrid   = 210
Wgrid   = 170

# Faille
Lfaille = 26.25
Wfaille = 21.25

# Enregistrer ou Visualiser
savefigures = False

# Nom du modele
modname = 'Niigata_Shibas_model'
# ********************************************************* #

# Nombre totale des points sources
NSR 	= Wgrid*Lgrid

# Lire les coordonnees de l'hypocentre (point de nucleation)
filename     = 'nucleationpoint.dat'
hypox, hypoy = readhypo(filename)

# Lire les coordonnes des points de grille et slip_max
filename   = 'slipdistribution.dat'
x, y, slip = readfileslip(NSR, filename)

# Tracer le slip
plotslip (Wgrid,Lgrid,Lfaille,Wfaille,slip,savefigures,modname+'_slip.png',hypox,hypoy)

# Lire les coordonnes des petites sources
filename     = 'subsources.dat'
xsub, ysub, radius = read_subsources (filename)

# Tracer les petites sources 
plotsubsources(xsub,ysub,radius,Lfaille,Wfaille,savefigures,modname+'_subsources.png',hypox,hypoy)

# Tracer la grille avec l'hypocentre
plotgrid(Wgrid,Lgrid,Lfaille,Wfaille,savefigures,modname+'_grid.png',hypox,hypoy,x,y)


print;
print '*** Done ***'
#