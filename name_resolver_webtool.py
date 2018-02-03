#Import the necessary packages
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np #Numpy maths
import pdb #Debugging

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
	
	#out = banyan_sigma_wrapper(name=name,ip='1.2.3',ra=data['RADEG'][0],dec=data['DECDEG'][0],pmra=data['PMRA'][0],pmdec=data['PMDEC'][0],epmra=data['ePMRA'][0],epmdec=data['ePMDEC'][0],rv=data['VRAD'][0],erv=data['eVRAD'][0],plx=data['PLX'][0],eplx=data['ePLX'][0])
	#pdb.set_trace()
	
	#Export to CSV file
	outfile = '/home/ipm/banyan/banyansigma/answer/info_'+name+'.dat'
	data[['Name','RADEG','DECDEG','PMRA','ePMRA','PMDEC','ePMDEC','REFPM','VRAD','eVRAD','PLX','ePLX']].to_csv(outfile,index=False)
