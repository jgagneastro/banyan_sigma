#Import the necessary packages
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np #Numpy maths
import pdb #Debugging
from banyan_sigma import banyan_sigma
from astropy.table import Table

def banyan_sigma_wrapper(name=None,ip=None,ra=None,dec=None,pmra=None,pmdec=None,epmra=None,epmdec=None,rv=None,erv=None,plx=None,eplx=None):
	
	#Launch the regular banyan_sigma
	stars_data = pd.DataFrame({'RA':ra,'DEC':dec,'PMRA':pmra,'EPMRA':epmra,'PMDEC':pmdec,'EPMDEC':epmdec,'RV':rv,'ERV':erv,'PLX':plx,'EPLX':plx},index=[0])
	use_rv = False
	use_plx = False
	if np.isfinite(stars_data['RV'][0]) and np.isfinite(stars_data['ERV'][0]):
		use_rv = True
	if np.isfinite(stars_data['PLX'][0]) and np.isfinite(stars_data['EPLX'][0]):
		use_plx = True
	output = banyan_sigma(Table.from_pandas(stars_data),use_rv=use_rv,use_plx=use_plx)
	#output = banyan_sigma(Table.from_pandas(stars_data),column_names={'RA':'RA','DEC':'DEC','PMRA':'PMRA','EPMRA':'EPMRA','PMDEC':'PMDEC','EPMDEC':'EPMDEC'},use_rv=use_rv,use_plx=use_plx)
	
	#Transform LN_P to 0-1 probabilities
	probs = np.exp(output['ALL'].values)
	output['ALL'].loc[0] = probs
	pdb.set_trace()
	
	#Save output probabilities to CSV
	outdir = '/home/gagne/www/banyansigma/answer/'
	output['ALL'].to_csv(outdir+'_prob'+name+'.dat')
	
	#Read all most probable RVs
	rv_opt = pd.DataFrame()
	erv_opt = pd.DataFrame()
	d_opt = pd.DataFrame()
	ed_opt = pd.DataFrame()
	for keys in output['ALL'].keys():
		rv_opt[keys] = [output[keys]['RV_OPT'][0]]
		erv_opt[keys] = [output[keys]['ERV_OPT'][0]]
		d_opt[keys] = [output[keys]['D_OPT'][0]]
		ed_opt[keys] = [output[keys]['ED_OPT'][0]]

	#Save optimal quantities to CSV files	
	d_opt.to_csv(outdir+'mdist_'+name+'.dat')
	ed_opt.to_csv(outdir+'emdist_'+name+'.dat')
	erv_opt.to_csv(outdir+'emvrad_'+name+'.dat')
	rv_opt.to_csv(outdir+'mvrad_'+name+'.dat')	

def name_resolver_webtool(name=None):
	
	if name is None:
		return
	
	#Add required columns
	Simbad.add_votable_fields('pmra','pmdec','pm_err_mina','pm_err_maja','pm_err_angle','pm_bibcode','plx','plx_error','rv_value','rvz_error')
	
	#Query Simbad
	result_table = Simbad.query_object(name)
	
	#Convert coordinates
	sky_pos = SkyCoord(result_table[0]['RA'],result_table[0]['DEC'], unit=(u.hour, u.deg), frame='icrs')
	
	#Convert proper motions
	epmra = np.sqrt((np.sin(result_table[0]['PM_ERR_ANGLE']*u.degree).value*result_table[0]['PM_ERR_MAJA'])**2 + (np.cos(result_table[0]['PM_ERR_ANGLE']*u.degree).value*result_table[0]['PM_ERR_MINA'])**2)
	epmdec = np.sqrt((np.cos(result_table[0]['PM_ERR_ANGLE']*u.degree).value*result_table[0]['PM_ERR_MAJA'])**2 + (np.sin(result_table[0]['PM_ERR_ANGLE']*u.degree).value*result_table[0]['PM_ERR_MINA'])**2)
	
	#Store in a pandas DataFrame for easy csv export
	data = pd.DataFrame({'Name':name,'RADEG':sky_pos.ra.degree,'DECDEG':sky_pos.dec.degree,'PMRA':result_table[0]['PMRA'],'ePMRA':epmra,'PMDEC':result_table[0]['PMDEC'],'ePMDEC':epmdec,'REFPM':result_table[0]['PM_BIBCODE'].decode('utf-8'),'VRAD':result_table[0]['RV_VALUE'],'eVRAD':result_table[0]['RVZ_ERROR'],'PLX':result_table[0]['PLX_VALUE'],'ePLX':result_table[0]['PLX_ERROR']},index=[0])
	
	#Make sure to pass NaNs when values are masked
	if (bool(result_table['RV_VALUE'].mask[0]) is True) or (bool(result_table['RVZ_ERROR'].mask[0]) is True):
		data['VRAD'] = np.nan
		data['eVRAD'] = np.nan
	if (bool(result_table['PLX_VALUE'].mask[0]) is True) or (bool(result_table['PLX_ERROR'].mask[0]) is True):
		data['PLX'] = np.nan
		data['ePLX'] = np.nan
	if (bool(result_table['PMRA'].mask[0]) is True) or (bool(result_table['PMDEC'].mask[0]) is True) or (bool(result_table['PM_ERR_MAJA'].mask[0]) is True) or (bool(result_table['PM_ERR_MINA'].mask[0]) is True) or (bool(result_table['PM_ERR_ANGLE'].mask[0]) is True):
		data['PMRA'] = np.nan
		data['ePMRA'] = np.nan
		data['PMDEC'] = np.nan
		data['ePMDEC'] = np.nan
	
	#TEST
	out = banyan_sigma_wrapper(name=name,ip='1.2.3',ra=data['RADEG'][0],dec=data['DECDEG'][0],pmra=data['PMRA'][0],pmdec=data['PMDEC'][0],epmra=data['ePMRA'][0],epmdec=data['ePMDEC'][0],rv=data['VRAD'][0],erv=data['eVRAD'][0],plx=data['PLX'][0],eplx=data['ePLX'][0])
	pdb.set_trace()
	
	#Export to CSV file
	outfile = '/home/gagne/www/banyansigma/answer/info_'+name+'.dat'
	data[['Name','RADEG','DECDEG','PMRA','ePMRA','PMDEC','ePMDEC','REFPM','VRAD','eVRAD','PLX','ePLX']].to_csv(outfile,index=False)

#import numpy as np #Numpy maths
#from scipy.special import erfc
#import os #Access to environment variables
#import pandas as pd #Pandas dataframes will be used to store BANYAN Sigma outputs
#from astropy.table import Table #Reading astro-formatted tables
#import warnings #Raise user-defined Python warnings
#import pdb #Debugging
#from scipy.stats import describe #Useful for debugging
#from scipy.misc import logsumexp #Useful to sum logarithms in a numerically stable way