__author__ = 'Jonathan Gagne'
__email__ = 'jonathan.gagne@astro.umontreal.ca'
__uri__ = "https://github.com/jgagneastro/banyan_sigma"
__license__ = "MIT"
__description__ = "A Bayesian classifier for young moving groups membership"

#Example code

#Import banyan_sigma
from banyan_sigma import *

#Define observables for an example star
ra=311.2911826481039
dec=-31.3425000799281

pmra=281.319
epmra=0.022
pmdec=-360.148
epmdec=0.019

plx=102.943
eplx=0.023

rv=-5.2
erv=0.7

#Determine membership probability
output = membership_probability(ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,epmra=epmra,epmdec=epmdec,plx=plx,eplx=eplx,rv=rv,erv=erv, use_plx=True, use_rv=True)

#Check the results
output.iloc[0]