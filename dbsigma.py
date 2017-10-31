from banyan_sigma import *
import numpy as np
from astropy.table import Table
import pdb
from scipy.stats import describe #Useful for debugging

#Debugging
def dbsigma():
	
	#Data file containing the parameters of Bayesian hypotheses
	parameters_file = os.path.dirname(__file__)+os.sep+'data'+os.sep+'banyan_sigma_parameters.fits'
	
	#Data file containing the performance metrics of young associations
	performance_file = os.path.dirname(__file__)+os.sep+'data'+os.sep+'banyan_sigma_metrics.fits'
	
	#Read the parameters of Bayesian hypotheses
	parameters_str = Table.read(parameters_file,format='fits')
	
	#Read the performance metrics of young associations
	
	#Read star data
	stars_str = Table.read('/Users/gagne/Downloads/DR7_tot_6_0.fits',format='fits')
	stars_str = stars_str[np.where((stars_str['PMRA'] != 0) & (stars_str['PMRA'] != -9999) & (stars_str['PMDEC'] != 0) & (stars_str['PMDEC'] != -9999))]
	
	column_data = {"EPMRA":"PMRAERR","EPMDEC":"PMDECERR"}
	test = banyan_sigma(stars_data=stars_str,column_names=column_data)
	
	test = banyan_sigma(ra=stars_str['RA'],dec=stars_str['DEC'],pmra=stars_str['PMRA'],pmdec=stars_str['PMDEC'],epmra=stars_str['PMRAERR'],epmdec=stars_str['PMDECERR'])
	pdb.set_trace()
	
	
	
	pdb.set_trace()
	test = banyan_sigma(stars_data=stars_str,column_data=column_data)
	
	pdb.set_trace()
	