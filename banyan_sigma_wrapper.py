#Import the necessary packages
import pandas as pd
import numpy as np #Numpy maths
from banyan_sigma import banyan_sigma
import pdb #Debugging

def banyan_sigma_wrapper(name=None,ip=None,ra=None,dec=None,pmra=None,pmdec=None,epmra=None,epmdec=None,rv=None,erv=None,plx=None,eplx=None):
	
	#Launch the regular banyan_sigma
	output = banyan_sigma(ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,epmra=epmra,epmdec=epmdec,rv=rv,erv=erv,plx=plx,eplx=eplx)
	
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