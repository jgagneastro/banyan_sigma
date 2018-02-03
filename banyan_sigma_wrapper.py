#Import the necessary packages
import pandas as pd
import numpy as np #Numpy maths
from banyan_sigma import banyan_sigma
import pdb #Debugging
#/home/ipm/banyan/python_launch_banyansigma.bash 'HIP9' '0.035320833333333336' '36.58595833333334' '-6.88' '0.5799999833106995' '8.57' '1.0399999618530273' '2007A&A...474..653V' '-9999' '-9999' '2.38' '0.9300000071525574'

def banyan_sigma_wrapper(name=None,ip=None,ra=None,dec=None,pmra=None,pmdec=None,epmra=None,epmdec=None,rv=None,erv=None,plx=None,eplx=None):
	
	#Parse missing values
	if rv == '-9999' or erv == '-9999':
		rv = np.nan
		erv = np.nan
	if plx == '-9999' or eplx == '-9999':
		plx = np.nan
		eplx = np.nan
	
	#Launch the regular banyan_sigma
	output = banyan_sigma(ra=float(ra),dec=float(dec),pmra=float(pmra),pmdec=float(pmdec),epmra=float(epmra),epmdec=float(epmdec),rv=float(rv),erv=float(erv),plx=float(plx),eplx=float(eplx))
	
	#Transform LN_P to 0-1 probabilities
	output_all = output['ALL']
	probs = output_all.iloc[0].values*1e2
	max_prob = 99.9
	probs = np.minimum(probs,max_prob)
	probs_formatted = np.round(np.minimum(probs,max_prob),1)
	output_all.loc[0] = probs_formatted
	print(output_all)
	
	#Save output probabilities to CSV
	outdir = '/home/ipm/banyan/banyansigma/answer/'
	output_all.to_csv(outdir+'prob_'+name+'.dat',index=False)
	
	#Read all most probable RVs
	rv_opt = pd.DataFrame()
	erv_opt = pd.DataFrame()
	d_opt = pd.DataFrame()
	ed_opt = pd.DataFrame()
	for keys in output['ALL'].keys():
		rv_opt[keys] = [np.round(output[keys]['RV_OPT'][0],1)]
		erv_opt[keys] = [np.round(output[keys]['ERV_OPT'][0],1)]
		d_opt[keys] = [np.round(output[keys]['D_OPT'][0],1)]
		ed_opt[keys] = [np.round(output[keys]['ED_OPT'][0],1)]
	
	print(rv_opt)
	#Save optimal quantities to CSV files	
	d_opt.to_csv(outdir+'mdist_'+name+'.dat',index=False)
	ed_opt.to_csv(outdir+'emdist_'+name+'.dat',index=False)
	erv_opt.to_csv(outdir+'emvrad_'+name+'.dat',index=False)
	rv_opt.to_csv(outdir+'mvrad_'+name+'.dat',index=False)