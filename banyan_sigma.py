#from scipy import constants as cs #constants of physics
#import emcee #MCMC
#import matplotlib.pyplot as plt #Plotting
#import math #Basic maths
#import numpy as np #More maths
#from astropy import constants as asc #Astronomical constants
#import datetime #Timing calculations
#import corner #Corner plots

import numpy as np #Numpy maths

#Maybe make this routine callable directly from the terminal
def banyan_sigma(period=2.425,eperiod=0.003,vsini=50.9,evsini=0.8,radii_bounds=[0.8,1.2],nMonte=int(1e4),nwalkers=12,nburnin=int(1e3),nThreads=4,figurename='corner.png'):
	#Period in hours
	#v sin i in km/s
	#Radii in R_Jup
	
	lnP = 0.
	return lnP
	
def banyan_sigma_solve_multivar(ra,dec,pmra,pdmec,epmra,epmdec,precision_matrix=None,center_vec=None,rv_measured=None,rv_error=None,dist_measured=None,dist_error=None,psira=None,psidec=None,epsira=None,epsidec=None,full_statistical_errors=False,lnP_only=False,kappa=None)
	"""Solve the radial velocity and distance marginalization integrals (if needed) and compute log(probability) with Bayes theroem for an array of stars and a single multivariate Gaussian XYZUVW model. This is a subroutine of banyan_sigma.
	
	multivar_model is IDL's "association_structure"
	#(ra,dec): Sky position (degrees)
	#(pmra,pmdec): Proper motion (mas/yr). pmra must include the cos(delta) term
	#(epmra,epmdec): Measurement errors on proper motion (mas/yr)
	#precision_matrix: Inverse of the covariance matrix [XYZUVW] of the multivariate Gaussian model (mixed units of pc and km/s)
	#center_vec: Central XYZUVW position of the multivariate Gaussian model (mixed units of pc and km/s)
	
	#(rv_measured,rv_error): Radial velocity measurement and error (km/s) - Optional inputs
	#(dist_measured,dist_error): Distance measurement and error (pc) - Optional inputs
	#(psira,psidec): Psi vector (described in Gagne et al., in preparation) describing the parallax motion of the star. This can be used to model the effect of parallax motion when a proper motion was measured from only two epochs ([mas/yr]) - Optional inputs
	#(epsira,epsidec): Measurement errors of the psi vector ([mas/yr]) - Optional inputs
	#full_statistical_errors: Compute [full statistical errors]
	#lnP_only: Only return the ln(probability)
	#kappa: [X]
	"""
	
	#Check for keyword consistency
	num_stars = np.size(ra)
	if np.size(dec) != num_stars or np.size(pmra) != num_stars or np.size(pmdec) != num_stars:
		void=1#Raise error message here
	
	