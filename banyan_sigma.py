"""
View the README.md file for a full description of this code and how to use it.
"""

#Import the necessary packages
import numpy as np #Numpy maths
from scipy.special import erfc
import os #Access to environment variables
import pandas as pd #Pandas dataframes will be used to store BANYAN Sigma outputs
from astropy.table import Table #Reading astro-formatted tables
import warnings #Raise user-defined Python warnings
import pdb #Debugging
from scipy.stats import describe #Useful for debugging
from scipy.special import logsumexp #Useful to sum logarithms in a numerically stable way

#A more user-friendly way to set break points
stop = pdb.set_trace

#A very small number used for numerical stability
tiny_number = 1e-318

#The total number of stars in the Besancon model within 300 pc to tranlate FPR to NFP
total_besancon_objects = 7152397.0

#Initiate some global constants
#1 AU/yr to km/s divided by 1000
kappa = 0.004743717361
#Not using "from astropy import units as u; kappa=u.au.to(u.km)/u.year.to(u.s)" because astropy defines one year as exactly 365.25 days instead of 365 days

#J2000.0 Equatorial position of the Galactic North (b=90 degrees) from Carrol and Ostlie
ra_pol = 192.8595
dec_pol = 27.12825

#J2000.0 Galactic latitude gb of the Celestial North pole (dec=90 degrees) from Carrol and Ostlie
l_north = 122.932

#Galactic Coordinates matrix
TGAL = (np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
	[0.4941094279, -0.4448296300, 0.7469822445],
	[-0.8676661490,  -0.1980763734, 0.4559837762]]))

#Initiate some secondary variables
sin_dec_pol = np.sin(np.radians(dec_pol))
cos_dec_pol = np.cos(np.radians(dec_pol))

#Main BANYAN_SIGMA routine
def banyan_sigma(stars_data=None,column_names=None,hypotheses=None,ln_priors=None,ntargets_max=1e7,ra=None,dec=None,pmra=None,pmdec=None,epmra=None,epmdec=None,dist=None,edist=None,rv=None,erv=None,psira=None,psidec=None,epsira=None,epsidec=None,plx=None,eplx=None,constraint_dist_per_hyp=None,constraint_edist_per_hyp=None,unit_priors=False,lnp_only=False,no_xyz=False,use_rv=None,use_dist=None,use_plx=None,use_psi=None,custom_models=None):
	
	#Automatically detect Astropy Tables and transform them to pandas dataframes
	if stars_data is not None:
		if isinstance(stars_data,Table):
			#First remove multi-dimensional columns to avoid crash
			for keys in stars_data.keys():
				if stars_data[keys].ndim != 1:
					stars_data.remove_column(keys)
			#Now transform to pandas dataframe
			stars_data = stars_data.to_pandas()
	
	#Check input consistency
	if stars_data is None and (ra is None or dec is None or pmra is None or pmdec is None or epmra is None or epmdec is None):
		raise ValueError('Either an input structure (stars_data) or all of the ra,dec,pmra,pmdec,epmra and epmdec keywords must be specified !')
	
	if constraint_dist_per_hyp is not None and constraint_edist_per_hyp is None:
		raise ValueError('f constraint_dist_per_hyp is specified, constraint_edist_per_hyp must also be specified !')
	
	#Default column names
	default_column_names = {'RA':'RA','DEC':'DEC','PMRA':'PMRA','PMDEC':'PMDEC','EPMRA':'EPMRA','EPMDEC':'EPMDEC'}
	if use_rv is True:
		default_column_names['RV'] = 'RV'
		default_column_names['ERV'] = 'ERV'
	if use_plx is True:
		default_column_names['PLX'] = 'PLX'
		default_column_names['EPLX'] = 'EPLX'
	if use_dist is True:
		default_column_names['DIST'] = 'DIST'
		default_column_names['EDIST'] = 'EDIST'
	if use_psi is True:
		default_column_names['PSIRA'] = 'PSIRA'
		default_column_names['PSIDEC'] = 'PSIDEC'
		default_column_names['EPSIRA'] = 'EPSIRA'
		default_column_names['EPSIDEC'] = 'EPSIDEC'
	
	#Merge user-issued column data with the default values (the user-issued values take predominance)
	if column_names is not None:
		column_names = {**default_column_names, **column_names}
	else:
		column_names = default_column_names
	
	#Check if a column named PLX, DIST, RV, PSIRA, etc. exist in stars_data but not in column_names. If this is the case, issue a warning so that the user understands that some data are not being considered.
	if stars_data is not None:
		if 'PLX' in stars_data.keys() and 'PLX' not in column_names.keys() and use_plx is None:
			warnings.warn('Parallaxes (PLX) were not read from the input data, because the PLX key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with the use_plx=True keyword to read them, or with use_plx=False to avoid this warning message.')
		if 'DIST' in stars_data.keys() and 'DIST' not in column_names.keys() and use_dist is None:
			warnings.warn('Distances (DIST) were not read from the input data, because the DIST key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with the use_dist=True keyword to read them, or with use_dist=False to avoid this warning message.')
		if 'RV' in stars_data.keys() and 'RV' not in column_names.keys() and use_rv is None:
			warnings.warn('Radial velocities (RV) were not read from the input data, because the RV key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with use_rv=True to read them, or with use_rv=False to avoid this warning message.')
		if ('PSIRA' in stars_data.keys() and 'PSIRA' not in column_names.keys()) or ('PSIDEC' in stars_data.keys() and 'PSIDEC' not in column_names.keys()) and use_psi is None:
			warnings.warn('The PSI parameters (PSIRA,PSIDEC) were not read from the input data, because the PSIRA and PSIDEC keys were not included in the column_data keyword of banyan_sigma(). You can also call banyan_sigma() with use_psi=True keyword to read them, or with use_psi=False to avoid this warning message.')
	
	#Create a table of data for BANYAN SIGMA to use
	if ra is not None:
		nobj = np.size(ra)
		zeros = np.zeros(nobj)
		data_table = pd.DataFrame({'RA':ra,'DEC':dec,'PMRA':pmra,'PMDEC':pmdec,'EPMRA':epmra,'EPMDEC':epmdec,'PSIRA':zeros,'PSIDEC':zeros,'EPSIRA':zeros,'EPSIDEC':zeros})
	if ra is None:
		nobj = stars_data.shape[0]
		zeros = np.zeros(nobj)
		data_table = pd.DataFrame({'RA':stars_data[column_names['RA']],'DEC':stars_data[column_names['DEC']],'PMRA':stars_data[column_names['PMRA']],'PMDEC':stars_data[column_names['PMDEC']],'EPMRA':stars_data[column_names['EPMRA']],'EPMDEC':stars_data[column_names['EPMDEC']],'PSIRA':zeros,'PSIDEC':zeros,'EPSIRA':zeros,'EPSIDEC':zeros})
	
	#Fill up the data table with stars_data if it is specified
	if stars_data is not None:
		for keys in column_names.keys():
			#Skip special keys
			if (keys == 'NAME') or (keys == 'PLX') or (keys == 'EPLX'):
				continue
			data_table[keys] = stars_data[column_names[keys]]
		if 'PLX' in column_names.keys():
			data_table['DIST'] = 1e3/stars_data[column_names['PLX']]
		if 'PLX' in column_names.keys() and 'EPLX' in column_names.keys():
			data_table['EDIST'] = 1e3/stars_data[column_names['PLX']]**2*stars_data[column_names['EPLX']]
	
	#Transform parallaxes to distances directly in data_table
	if 'PLX' in data_table.keys() and 'EPLX' in data_table.keys():
		data_table['EDIST'] = 1e3/data_table['PLX']**2*data_table['EPLX']
		data_table = data_table.drop('EPLX', 1)
	if 'PLX' in data_table.keys():
		data_table['DIST'] = 1e3/data_table['PLX']
		data_table = data_table.drop('PLX', 1)
	
	#If measurements are specified as keywords, put them in the data table
	if ra is not None:
		data_table['RA'] = ra
	if dec is not None:
		data_table['DEC'] = dec
	if pmra is not None:
		data_table['PMRA'] = pmra
	if pmdec is not None:
		data_table['PMDEC'] = pmdec
	if epmra is not None:
		data_table['EPMRA'] = epmra
	if epmdec is not None:
		data_table['EPMDEC'] = epmdec
	if plx is not None:
		data_table['DIST'] = 1e3/plx
	if plx is not None and eplx is not None:
		data_table['EDIST'] = 1e3/plx**2*eplx
	if dist is not None:
		data_table['DIST'] = dist
	if edist is not None:
		data_table['EDIST'] = edist
	if rv is not None:
		data_table['RV'] = rv
	if erv is not None:
		data_table['ERV'] = erv
	if psira is not None:
		data_table['PSIRA'] = psira
	if psidec is not None:
		data_table['PSIDEC'] = psidec
	if epsira is not None:
		data_table['EPSIRA'] = epsira
	if epsidec is not None:
		data_table['EPSIDEC'] = epsidec
	
	#Check for unphysical data
	if np.max((data_table['RA'] < 0.) | (data_table['RA'] >= 360.)) != 0:
		raise ValueError('Some RA values are unphysical')
	if np.max((data_table['DEC'] < -90.) | (data_table['DEC'] > 90.)) != 0:
		raise ValueError('Some DEC values are unphysical')
	if np.max((data_table['EPMRA'] < 0.) | (data_table['EPMDEC'] < 0.)) != 0:
		raise ValueError('Some EPMRA or EPMDEC values are unphysical')
	if np.max((np.isnan(data_table['RA']) | (np.isnan(data_table['DEC'])) | (np.isnan(data_table['PMRA'])) | (np.isnan(data_table['PMDEC'])) | (np.isnan(data_table['EPMRA'])) | (np.isnan(data_table['EPMDEC'])))) != 0:
		raise ValueError('The observables ra,dec,pmra,pmdec,epmra and epmdec must be specified (and finite) for each object !')
	if 'RV' in data_table.keys() and 'ERV' not in data_table.keys():
		raise ValueError('RV is defined in the data table but not ERV')
	if 'DIST' in data_table.keys() and 'EDIST' not in data_table.keys():
		raise ValueError('DIST is defined in the data table but not EDIST')
	if 'ERV' in data_table.keys():
		if np.max(data_table['ERV'] <= 0.):
			raise ValueError('Some ERV values are unphysical')
	if 'RV' in data_table.keys() and 'ERV' in data_table.keys():
		if np.max(np.isfinite(data_table['RV']) & np.isnan(data_table['ERV'])):
			raise ValueError('Some RV values are specified without ERV')
	if 'DIST' in data_table.keys() and 'EDIST' in data_table.keys():
		if np.max((data_table['DIST'] < 0.) | (data_table['EDIST'] <= 0.)):
			raise ValueError('Some DIST or EDIST values are unphysical')
		if np.max(np.isfinite(data_table['DIST']) & np.isnan(data_table['EDIST'])):
			raise ValueError('Some DIST values are specified without EDIST')
	if np.max(((data_table['PSIRA'] != 0.) | (data_table['PSIDEC'] != 0.)) & ((data_table['EPSIRA'] == 0.) | (data_table['EPSIDEC'] == 0.)) | (data_table['EPSIRA'] < 0.) | (data_table['EPSIDEC'] < 0.)):
			raise ValueError('Some EPSIRA or EPSIDEC values are unphysical')
	
	#Fill the data table with empty RVs and distances if they were not specified
	if 'RV' not in data_table.keys():
		data_table['RV'] = np.nan
	if 'ERV' not in data_table.keys():
		data_table['ERV'] = np.nan
	if 'DIST' not in data_table.keys():
		data_table['DIST'] = np.nan
	if 'EDIST' not in data_table.keys():
		data_table['EDIST'] = np.nan
	
	if custom_models is not None:
		parameters_str = custom_models
	else:
		#Data file containing the parameters of Bayesian hypotheses
		parameters_file = os.path.dirname(__file__)+os.sep+'data'+os.sep+'banyan_sigma_parameters.fits'
		
		#Check if the file exists
		if not os.path.isfile(parameters_file):
			raise ValueError('The multivariate Gaussian parameters file could not be found ! Please make sure that you did not move "'+os.sep+'data'+os.sep+'banyan_sigma_parameters.fits" from the same path as the Python file banyan_sigma.py !')
		
		#Read the parameters of Bayesian hypotheses
		parameters_str = Table.read(parameters_file,format='fits')
	
	#Remove white spaces in names
	parameters_str['NAME'] = np.chararray.strip(np.array(parameters_str['NAME']))

	#Index the table by hypothesis name
	parameters_str.add_index('NAME')
	npar = np.size(parameters_str)
	
	#Build a unique list of Bayesian hypotheses
	if hypotheses is None:
		hypotheses = np.array(parameters_str['NAME'])
		indexes = np.unique(hypotheses,return_index=True)[1]
		hypotheses = hypotheses[sorted(indexes)]
	
	#Make sure that hypotheses are all upper case
	#Also make sure that all hypotheses are not in bytes format
	hypotheses = np.array([hyp.upper().decode('UTF-8') for hyp in hypotheses.tolist()])
	nhyp = hypotheses.size
	
	#If constraint_dist_per_hyp is set, check that all hypotheses are included
	if constraint_dist_per_hyp is not None:
		if sorted(constraint_dist_per_hyp.keys()) != sorted(constraint_edist_per_hyp.keys()):
			raise ValueError('The tag names of constraint_dist_per_hyp and constraint_edist_per_hyp are different')
		if sorted(constraint_dist_per_hyp.keys()) != sorted(hypotheses.tolist()):
			raise ValueError('The tag names of constraint_dist_per_hyp and the list of Bayesian hypotheses are different')
		
		#Build constraint_dist_per_hyp into an array
		dist_per_hyp_arr = np.empty((nobj,nhyp))*np.nan
		edist_per_hyp_arr = np.empty((nobj,nhyp))*np.nan
		#Read the distance constraints for each Bayesian hypothesis
		for i in range(nhyp):
			dist_per_hyp_arr[:,i] = constraint_dist_per_hyp[hypotheses[i]]
			edist_per_hyp_arr[:,i] = constraint_edist_per_hyp[hypotheses[i]]
		
		#Verify that all distance constraints are physical
		if np.max(dist_per_hyp_arr < 0. | edist_per_hyp_arr <= 0.):
			raise ValueError('Some of the specified constraint_dist_per_hyp or constraint_edist_per_hyp values are unphysical')
		if np.max(np.isfinite(dist_per_hyp_arr) & np.isnan(edist_per_hyp_arr)):
			raise ValueError('Some of the specified constraint_edist_per_hyp are not finite where constraint_dist_per_hyp are finite')
		
		#Check that either all or none of the distance constraints are finite for a given object
		if np.max(np.isfinite(np.nansum(dist_per_hyp_arr,axis=1)) and np.isnan(np.sum(dist_per_hyp_arr,axis=1))):
			raise ValueError('The constraint_dist_per_hyp and constraint_edist_per_hyp values must be all finite or all non-finite for a given star')
	
	#Override priors to unity if the keyword unit_priors is set
	if unit_priors is True:
		parameters_str['LN_PRIOR'] = 0.
	
	#Determine whether a trigonometric distance or a per-hypothesis distance constraint was set
	if constraint_dist_per_hyp is not None:
		distance_is_set = (np.isfinite(data_table['DIST']) | np.isfinite(np.nansum(dist_per_hyp_arr,axis=1)))
	else:
		distance_is_set = np.isfinite(data_table['DIST'])
	
	#Assign the correct Bayesian priors to each star
	g_pm = (np.where(np.isnan(data_table['RV']) & (~distance_is_set)))[0]
	g_pm_rv = (np.where(np.isfinite(data_table['RV']) & (~distance_is_set)))[0]
	g_pm_dist = (np.where(np.isnan(data_table['RV']) & distance_is_set))[0]
	g_pm_rv_dist = (np.where(np.isfinite(data_table['RV']) & distance_is_set))[0]
	ln_priors_nd = np.zeros((nobj,nhyp))
	ln_priors_nd_manual = np.zeros((nobj,nhyp))
	for i in range(nhyp):
		#Skip the field hypotheses as they do not have a Bayesian prior
		if hypotheses[i].find('FIELD') != -1:
			continue
		#Read the parameters structure to identify the 4 priors associated with a given young association
		ln_priors_i = parameters_str.loc[hypotheses[i]]['LN_PRIOR']
		#In the cases where only one prior is designated, assign it to all stars
		if ln_priors_i.size == 1:
			ln_priors_nd[:,i] = ln_priors_i[0]
		else:
			#Otherwise assign them properly as a function of available observables
			ln_priors_nd[g_pm,i] = ln_priors_i[0]
			ln_priors_nd[g_pm_rv,i] = ln_priors_i[1]
			ln_priors_nd[g_pm_dist,i] = ln_priors_i[2]
			ln_priors_nd[g_pm_rv_dist,i] = ln_priors_i[3]
	
	#Include manual priors if they are specified as an input structure
	if ln_priors is not None:
		for i in range(nhyp):
			#The field hypotheses *can* have manual priors
			if hypotheses[i] not in ln_priors.keys():
				warnings.warn('The prior for hypothesis '+hypotheses[i]+' was left to its default value as it was not specified manually')
				continue
			ln_priors_nd_manual[:,i] = ln_priors[hypotheses[i]]
		
		#Normalize manual priors with the field hypothesis (because they get applied only on young associations)
		gnorm = np.where(['FIELD' in hyp for hyp in hypotheses.tolist()])
		norm_priors_1d = logsumexp(ln_priors_nd_manual[:,gnorm[0]],axis=1)
		ln_priors_nd_manual -= np.tile(norm_priors_1d,(nhyp,1)).transpose()
		
		#Apply the manual priors on top of the default priors
		ln_priors_nd += ln_priors_nd_manual
	
	#If both trigonometric distances and per-hypothesis distance constraints are set, transform the per-hypothesis distance constraints into priors
	both_distances_set = []
	if constraint_dist_per_hyp is not None:
		both_distances_set = np.where(np.isfinite(data_table['DIST']) & np.isfinite(np.nansum(dist_per_hyp_arr,axis=1)))
	if np.size(both_distances_set) != 0:
		xdist_measured = np.tile(data_table['DIST'].iloc[both_distances_set[0]],(nhyp,1)).transpose()
		xedist_measured = np.tile(data_table['EDIST'].iloc[both_distances_set[0]],(nhyp,1)).transpose()
		ln_prob_dist_differences = -(xdist_measured-dist_per_hyp_arr[both_distances_set[0],:])**2/(2.0*(xedist_measured**2+edist_per_hyp_arr[both_distances_set[0],:]**2))
		
		#Treat these values as priors so normalize them with the field hypotheses (because they get applied only on young associations)
		gnorm = np.where(['FIELD' in hyp for hyp in hypotheses.tolist()])
		norm_priors_1d = logsumexp(ln_prob_dist_differences[:,gnorm[0]],axis=1)
		ln_prob_dist_differences -= np.tile(norm_priors_1d,(nhyp,1)).transpose()
		
		#Apply these values on the priors
		ln_priors_nd[both_distances_set[0],L] += ln_prob_dist_differences
		
		#Remove the per-hypothesis distance constraints on these particular objects and just keep the trigonometric distances
		dist_per_hyp_arr[both_distances_set[0],:] = np.nan
		edist_per_hyp_arr[both_distances_set[0],:] = np.nan
	
	#Initiate an array that will contain the ln probabilities if those are the only required outputs
	if lnp_only is True:
		all_lnprobs = np.empty((nobj,nhyp))*np.nan
	
	#Loop on hypotheses to run BANYAN Sigma on
	output_str_allhyps_list = []
	for i in range(nhyp):
		#print("HYP "+str(i))
		
		#If constraint_dist_per_hyp is set, determine which distance constraint must be used now
		dist_for_this_hypothesis = data_table['DIST'].values
		edist_for_this_hypothesis = data_table['EDIST'].values
		if constraint_dist_per_hyp is not None:
			gdist_per_hyp = np.where(np.isfinite(dist_per_hyp_arr[:,i]))
			dist_for_this_hypotheses[gdist_per_hyp[0]] = dist_per_hyp[gdist_per_hyp[0],i]
			edist_for_this_hypotheses[gdist_per_hyp[0]] = edist_per_hyp_arr[gdist_per_hyp[0],i]
		
		#Loop over individual multivariate Gaussians if the model is a mixture
		ngauss = np.size(parameters_str.loc[hypotheses[i]])
		
		output_str_multimodel_list = []
		if lnp_only is True:
			all_lnprobs_hypi = np.zeros((nobj,ngauss))
		
		for gaussi in range(ngauss):
			
			#Somehow we cannot access the Gaussian index without the table breaking when there is just one Gaussian component, so here we grab the right table row
			if ngauss == 1:
				parameters_str_row = parameters_str.loc[hypotheses[i]]
			else:
				parameters_str_row = parameters_str.loc[hypotheses[i]][gaussi]
			
			#Determine how many batches will be needed to avoid saturating the RAM
			nbatches = np.int(np.ceil(nobj/ntargets_max))
			output_str_list = []
			for ci in range(nbatches):
				#Determine the indices of the stars to be selected
				ind_from = np.int(np.round(ci*ntargets_max))
				ind_to = np.int(ind_from + np.round(ntargets_max))
				ind_to = np.minimum(ind_to,np.int(nobj))
				
				#Create a sub-structure of input data
				data_table_ci = data_table[ind_from:ind_to]
				dist_for_this_hypothesis_ci = dist_for_this_hypothesis[ind_from:ind_to]
				edist_for_this_hypothesis_ci = edist_for_this_hypothesis[ind_from:ind_to]
				nobj_ci = np.size(data_table_ci)
				
				#Solve the BANYAN Sigma integrals for this hypothesis and this batch of targets
				output_str_ci = banyan_sigma_solve_multivar(data_table_ci['RA'].values,data_table_ci['DEC'].values,data_table_ci['PMRA'].values,data_table_ci['PMDEC'].values,data_table_ci['EPMRA'].values,data_table_ci['EPMDEC'].values,rv_measured=data_table_ci['RV'].values,rv_error=data_table_ci['ERV'].values,dist_measured=dist_for_this_hypothesis_ci,dist_error=edist_for_this_hypothesis_ci,psira=data_table_ci['PSIRA'].values,psidec=data_table_ci['PSIDEC'].values,psira_error=data_table_ci['EPSIRA'].values,psidec_error=data_table_ci['EPSIDEC'].values,precision_matrix=parameters_str_row['PRECISION_MATRIX'],center_vec=parameters_str_row['CENTER_VEC'],precision_matrix_determinant=parameters_str_row['PRECISION_DETERM'])
				
				#Store the log of probabilities if those are the only required output
				if lnp_only is True:
					all_lnprobs_hypi[ind_from:ind_to,gaussi] = output_str_ci['LN_P']
					continue
				
				#Append the dataframe in the Python list
				output_str_list.append(output_str_ci)
			
			#Contatenate the list of Dataframes
			output_str = pd.concat(output_str_list,ignore_index=True)
			
			#Reformat the output structure if this hypothesis is a multivariate Gaussian mixture
			if ngauss != 1:
				#Use column multi-indexing to add a second title to the columns, which corresponds to the ID if the Gaussian mixture component
				dataframe_column_names = output_str.columns
				output_str.columns = [np.array(dataframe_column_names),np.array(np.tile('Gauss'+str(gaussi),dataframe_column_names.size))]
				output_str_multimodel_list.append(output_str)
		
		#If only log probs are required, compile them in the main array
		if lnp_only is True:
			if ngauss == 1:
				all_lnprobs[:,i] = all_lnprobs_hypi
			else:
				weights = parameters_str.loc[hypotheses[i]]['COEFFICIENT']
				weights /= np.sum(weights)
				all_lnprobs[:,i] = logsumexp(np.tile(np.log(weights),(nobj,1))+all_lnprobs_hypi,axis=1)
			continue
		
		#Reformat the output structure if there is more than one multivariate gaussian
		if ngauss != 1:
			#Concatenate the list of pandas dataframes into a single dataframe
			output_str_multimodel = pd.concat(output_str_multimodel_list,axis=1)
			
			#Create a 2D array of weights to combine the Gaussian mixture components
			weights = parameters_str.loc[hypotheses[i]]['COEFFICIENT']
			weights /= np.sum(weights)
			logweights_2d = np.tile(np.log(weights),(nobj,1))
			
			#Combine each column of the dataframe with a weighted average
			output_str = pd.DataFrame()
			#Had to add a .values here
			for coli in output_str_multimodel.columns.get_level_values(0):
				output_str[coli] = logsumexp(logweights_2d+output_str_multimodel[coli].values,axis=1)
	
		#Use column multi-indexing to add a second title to the columns, which corresponds to the name of the Bayesian hypothesis
		dataframe_column_names = output_str.columns
		output_str.columns = [np.array(dataframe_column_names),np.array(np.tile(hypotheses[i],dataframe_column_names.size))]
		
		#Add the dataframe to the per-hypothesis list of dataframes
		output_str_allhyps_list.append(output_str)
	
	#Concatenate the list of pandas dataframes into a single dataframe
	output_str_all = pd.concat(output_str_allhyps_list,axis=1)
	
	#Fetch all log probabilities (if lnp_only is set, this variable already exists)
	if lnp_only is False:
		all_lnprobs = output_str_all['LN_P'].values
	
	#Normalize probabilities directly in log space
	ln_norm_output = all_lnprobs - np.tile(logsumexp(all_lnprobs,axis=1),(nhyp,1)).transpose()
	
	#Compute [0,1] probabilities
	norm_output = np.exp(ln_norm_output)
	
	#Identify hypotheses that correspond to moving groups or associations
	yind = (np.where(np.array([hypothesis.find('FIELD') == -1 for hypothesis in hypotheses])))[0]
	
	#Create an array of normalized YMG probabilities (no field)
	ln_norm_output_only_ymg = all_lnprobs[:,yind] - np.tile(logsumexp(all_lnprobs[:,yind],axis=1),(yind.size,1)).transpose()
	
	#Calculate the weighted YMG prior
	ln_prior_moving_groups = logsumexp(ln_priors_nd[:,yind]+ln_norm_output_only_ymg,axis=1)
	
	#Identify hypotheses that correspond to the field
	ffind = (np.where(np.array([hypothesis.find('FIELD') != -1 for hypothesis in hypotheses])))[0]
	
	#Weight the priors w/r/t the Bayesian probabilities and project these priors onto the field. This is a way to avoid having the priors change the relative moving group probabilities, as their goal is strictly to maximize young association vs FIELD classification performance
  	#Normalize probabilities directly in log space, projecting the inverse young association prior on the field probability
	ln_P_with_prior = all_lnprobs
	ln_P_with_prior[:,ffind] -= np.tile(ln_prior_moving_groups,(ffind.size,1)).transpose()
	#Renormalize
	ln_norm_output_prior = ln_P_with_prior - np.tile(logsumexp(ln_P_with_prior,axis=1),(nhyp,1)).transpose()
	
	#Return log probabilities if this is the only required output
	if lnp_only is True:
		return ln_norm_output_prior
	
	#Compute [0,1] probabilities
	norm_output_prior = np.exp(ln_norm_output_prior)
	
	#Data file containing the parameters of Bayesian hypotheses
	metrics_computed = False
	metrics_file = os.path.dirname(__file__)+os.sep+'data'+os.sep+'banyan_sigma_metrics.fits'
	
	#Check if the file exists
	if not os.path.isfile(metrics_file):
		warnings.warn('The performance metrics file could not be found ! Performance metrics will not be calculated. Please make sure that you did not move "'+os.sep+'data'+os.sep+'banyan_sigma_metrics.fits" from the same path as the Python file banyan_sigma.py !')
	
	#Avoid computing biased metrics if the unit_priors keyword was set
	if os.path.isfile(metrics_file) and unit_priors is False:
		metrics_str = Table.read(metrics_file,format='fits')
		#Remove white spaces in association names
		metrics_str['NAME'] = np.chararray.strip(np.array(metrics_str['NAME']))
		#Index the table by hypothesis name
		metrics_str.add_index('NAME')
		
		#Loop on young associations to determine their individual metrics
		tpr = np.empty((nobj,yind.size))*np.nan
		fpr = np.empty((nobj,yind.size))*np.nan
		ppv = np.empty((nobj,yind.size))*np.nan
		for yindi in range(yind.size):
			#Calculate the individual normalized probabilities for a given young association
			probs_yindi = np.exp(ln_norm_output_prior[:,yindi] - logsumexp(ln_norm_output_prior[:,[yindi,ffind[0]]],axis=1))
			#Store the interpolated values depending on observables
			if g_pm.size != 0:
				mode_index = 0
				tpr[g_pm,yindi] = np.interp(probs_yindi[g_pm],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['TPR'][mode_index,:])
				fpr[g_pm,yindi] = np.interp(probs_yindi[g_pm],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['FPR'][mode_index,:])
				ppv[g_pm,yindi] = np.interp(probs_yindi[g_pm],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['PPV'][mode_index,:])
			if g_pm_rv.size != 0:
				mode_index = 1
				tpr[g_pm_rv,yindi] = np.interp(probs_yindi[g_pm_rv],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['TPR'][mode_index,:])
				fpr[g_pm_rv,yindi] = np.interp(probs_yindi[g_pm_rv],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['FPR'][mode_index,:])
				ppv[g_pm_rv,yindi] = np.interp(probs_yindi[g_pm_rv],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['PPV'][mode_index,:])
			if g_pm_dist.size != 0:
				mode_index = 2
				tpr[g_pm_dist,yindi] = np.interp(probs_yindi[g_pm_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['TPR'][mode_index,:])
				fpr[g_pm_dist,yindi] = np.interp(probs_yindi[g_pm_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['FPR'][mode_index,:])
				ppv[g_pm_dist,yindi] = np.interp(probs_yindi[g_pm_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['PPV'][mode_index,:])
			if g_pm_rv_dist.size != 0:
				mode_index = 3
				tpr[g_pm_rv_dist,yindi] = np.interp(probs_yindi[g_pm_rv_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['TPR'][mode_index,:])
				fpr[g_pm_rv_dist,yindi] = np.interp(probs_yindi[g_pm_rv_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['FPR'][mode_index,:])
				ppv[g_pm_rv_dist,yindi] = np.interp(probs_yindi[g_pm_rv_dist],metrics_str.loc[hypotheses[yind[yindi]]]['PROBS'],metrics_str.loc[hypotheses[yind[yindi]]]['PPV'][mode_index,:])
		
		#Build the combination weights
		ln_weights = np.copy(ln_norm_output_only_ymg)
		#Any group with less than 1% probability is ignored to avoid propagating potential NaNs
		ln_weights[np.where(ln_weights < np.log(1e-2))] = np.log(tiny_number)
		#Re-normalize weights
		ln_weights -= np.tile(logsumexp(ln_weights,axis=1),(yind.size,1)).transpose()
		
		#Calculate the weighted metrics
		tpr_weighted = np.exp(logsumexp(np.log(np.maximum(tpr,tiny_number))+ln_weights,axis=1))
		fpr_weighted = np.exp(logsumexp(np.log(np.maximum(fpr,tiny_number))+ln_weights,axis=1))
		ppv_weighted = np.exp(logsumexp(np.log(np.maximum(ppv,tiny_number))+ln_weights,axis=1))
		metrics_computed = True
	
	#Determine the most probable hypothesis
	most_probable_index = np.nanargmax(norm_output_prior,axis=1)
	
	#Loop on objects to determine lists of good hypotheses
	hyp_lists = []
	best_ya = []
	norm_output_only_ymg = np.exp(ln_norm_output_only_ymg)
	for obji in range(nobj):
		#Identify all young associations with relative P>5%
		ind_obji = (np.where(norm_output_only_ymg[obji,:] > .05))[0]
		if len(ind_obji) == 0:
			hyp_lists.append('FIELD')
			best_ya.append('FIELD')
			continue
		
		#Find the most probable moving group
		best_ya_ind = np.nanargmax(norm_output_only_ymg[obji,:])
		best_ya.append(hypotheses[yind][best_ya_ind])
		
		#Sort by decreasing P
		ind_obji = ind_obji[np.flip(np.argsort(norm_output_only_ymg[obji,ind_obji]),axis=0)]
		#Build a list of associations
		if len(ind_obji) > 1:
			hyp_lists.append(';'.join([x+y for x,y in zip(hypotheses[yind][ind_obji].tolist(),['('+str(x)+')' for x in np.round(norm_output_only_ymg[obji,ind_obji]*1e2).astype(int).tolist()])]))
		if len(ind_obji) == 1:
			hyp_lists.append(hypotheses[yind][best_ya_ind])
	
	#Build a final output dataframe
	output_final = pd.DataFrame()
	
	#Store the star names if they are given
	if 'NAME' in data_table.keys():
		output_final['NAME'] = data_table['NAME']
	
	#Store global results
	output_final['YA_PROB'] = np.nansum(norm_output_prior[:,yind],axis=1)
	output_final['LIST_PROB_YAS'] = hyp_lists
	output_final['BEST_HYP'] = hypotheses[most_probable_index]
	output_final['BEST_YA'] = best_ya
	
	#Add a second column title "General"
	dataframe_column_names = output_final.columns
	output_final.columns = [np.array(dataframe_column_names),np.array(np.tile('Global',dataframe_column_names.size))]
	
	if metrics_computed is True:
		output_final['TPR','Metrics'] = tpr_weighted
		output_final['FPR','Metrics'] = fpr_weighted
		output_final['PPV','Metrics'] = ppv_weighted
		output_final['NFP','Metrics'] = fpr_weighted*total_besancon_objects
	
	#Create a Dataframe with all probabilities
	probs_frame = pd.DataFrame(norm_output_prior,columns=[np.array(np.tile('ALL',nhyp)),hypotheses])
	
	#Add the per-group stuff
	if metrics_computed is True:
		output_final = pd.concat([output_str_all.swaplevel(axis=1),probs_frame,output_final.swaplevel(axis=1)[['Metrics']],output_final.swaplevel(axis=1)[['Global']].swaplevel(axis=1)],axis=1)
	else:
		output_final = pd.concat([output_str_all.swaplevel(axis=1),probs_frame,output_final.swaplevel(axis=1)[['Global']].swaplevel(axis=1)],axis=1)
	
	#Add star names if they were provided
	if 'NAME' in data_table.keys():
		output_final.index = data_table['NAME']
	
	#Return the final structure
	return output_final
	
def banyan_sigma_solve_multivar(ra,dec,pmra,pmdec,pmra_error,pmdec_error,precision_matrix=None,center_vec=None,rv_measured=None,rv_error=None,dist_measured=None,dist_error=None,psira=None,psidec=None,psira_error=None,psidec_error=None,lnP_only=False,precision_matrix_determinant=None,debug=False):
	#PROBLEM: PSIRA_ERROR AND PSIDEC_ERROR ARE NOT USED ?
	"""
	Solve the radial velocity and distance marginalization integrals (if needed) and compute log(probability) with Bayes theorem for an array of stars and a single multivariate Gaussian XYZUVW model. This is a subroutine of banyan_sigma.
	
	Temporary note: multivar_model is IDL's "association_structure"
	
	params (ra,dec): Sky position (degrees)
	params (pmra,pmdec): Proper motion (mas/yr). pmra must include the cos(delta) term
	params (pmra_error,pmdec_error): Measurement errors on proper motion (mas/yr)
	param precision_matrix: Inverse of the covariance matrix [XYZUVW] of the multivariate Gaussian model (mixed units of pc and km/s)
	param precision_matrix_determinant; [X]
	param center_vec: Central XYZUVW position of the multivariate Gaussian model (mixed units of pc and km/s)
	params (rv_measured,rv_error): Radial velocity measurement and error (km/s) - Optional inputs
	params (dist_measured,dist_error): Distance measurement and error (pc) - Optional inputs
	params (psira,psidec): Psi vector (described in Gagne et al., in preparation) describing the parallax motion of the star. This can be used to model the effect of parallax motion when a proper motion was measured from only two epochs ([mas/yr]) - Optional inputs
	params (epsira,epsidec): Measurement errors of the psi vector ([mas/yr]) - Optional inputs
	keyword full_statistical_errors: Compute [full statistical errors]
	keyword lnP_only: Only return the ln(probability)
	"""
	
	#Check for parameter consistency
	num_stars = np.size(ra)
	if np.size(dec) != num_stars or np.size(pmra) != num_stars or np.size(pmdec) != num_stars or np.size(pmra_error) != num_stars or np.size(pmdec_error) != num_stars:
		raise ValueError('The dimensions ra, dec, pmra, pmdec, pmra_error and pmdec_error do not agree. They must all be numpy arrays of the same length.')
	
	#Check for radial velocity keyword consistencies
	if rv_measured is not None or rv_error is not None:
		if np.size(rv_measured) != num_stars or np.size(rv_error) != num_stars:
			raise ValueError('The dimensions of rv_measured or rv_error do not agree with those of ra, etc. They must all be numpy arrays of the same length.')
	
	#Check for distance keyword consistencies
	if dist_measured is not None or dist_error is not None:
		if np.size(dist_measured) != num_stars or np.size(dist_error) != num_stars:
			raise ValueError('The dimensions of dist_measured or dist_error do not agree with those of ra, etc. They must all be numpy arrays of the same length.')
	
	#Check for psi keyword consistencies
	if psira is not None or psidec is not None or psira_error is not None or psidec_error is not None:
		if np.size(psira) != num_stars or np.size(psidec) != num_stars or np.size(psira_error) != num_stars or np.size(psidec_error) != num_stars:
			raise ValueError('The dimensions of psira, psidec, psira_error or psidec_error do not agree with those of ra, etc. They must all be numpy arrays of the same length.')
	
	#Check that center_vec is a 6-elements array
	if np.shape(center_vec) != (6,):
		raise ValueError('center_vec must be a 6-elements numpy array.')
	
	#Check that precision_matrix is a 6x6 matrix
	if np.shape(precision_matrix) != (6, 6):
		raise ValueError('precision_matrix must be a 6x6-elements numpy array.')
	
	#Compute Galactic coordinates
	(gl,gb) = equatorial_galactic(ra,dec)
	
	#lambda is defined in Gagne et al. (2017, ApJS, X, Y, equation 7)
	cos_gl = np.cos(np.radians(gl))
	cos_gb = np.cos(np.radians(gb))
	sin_gl = np.sin(np.radians(gl))
	sin_gb = np.sin(np.radians(gb))
	lambda_vector = np.array([cos_gb*cos_gl,cos_gb*sin_gl,sin_gb]).transpose()
	
	#Build matrices A and B to convert sky quantities in the Galactic coordinates frame. The A matrix is defined in Gagne et al. (2017, ApJS, X, Y, equation 7)
	A_matrix = np.zeros((num_stars,3,3))
	cos_ra = np.cos(np.radians(ra))
	cos_dec = np.cos(np.radians(dec))
	sin_ra = np.sin(np.radians(ra))
	sin_dec = np.sin(np.radians(dec))
	A_matrix[:,0,0] = cos_ra * cos_dec
	A_matrix[:,1,0] = sin_ra * cos_dec
	A_matrix[:,2,0] = sin_dec
	A_matrix[:,0,1] = -sin_ra
	A_matrix[:,1,1] = cos_ra
	A_matrix[:,0,2] = -cos_ra * sin_dec
	A_matrix[:,1,2] = -sin_ra * sin_dec
	A_matrix[:,2,2] = cos_dec
	
	#The B matrix is not directly referenced in the BANYAN Sigma paper.
	B_matrix = matrix_set_product_A_single(TGAL,A_matrix)
	
	#The M vector is defined in Gagne et al. (2017, ApJS, X, Y, equation 7)
	M_vector = matrix_vector_set_product_v_single(B_matrix,np.array([1.0,0.0,0.0]))
	
	#The N vector is defined in Gagne et al. (2017, ApJS, X, Y, equation 7)
	N_vector_sub = np.array([np.zeros(num_stars), np.array(kappa*pmra), np.array(kappa*pmdec)]).transpose()
	N_vector = matrix_vector_set_product(B_matrix,N_vector_sub)
	
	#The varphi vector is defined in Gagne et al. (2017, ApJS, X, Y, equation 20)
	if psira is not None:
		varphi_vector_sub = np.array([np.zeros(num_stars),np.array(kappa*psira), np.array(kappa*psidec)]).transpose()
		varphi_vector = matrix_vector_set_product(B_matrix,varphi_vector_sub)
	
	#OMEGA is defined in Gagne et al. (2017, ApJS, X, Y, equation 6)
	zero_vector = np.zeros([num_stars,3])
	OMEGA_vector = np.concatenate((zero_vector,M_vector),axis=1)
	
	#GAMMA is defined in Gagne et al. (2017, ApJS, X, Y, equation 6)
	GAMMA_vector = np.concatenate((lambda_vector,N_vector),axis=1)
	
	#PHI is defined in Gagne et al. (2017, ApJS, X, Y, equation 20)
	if psira is not None:
		PHI_vector = np.concatenate((zero_vector,varphi_vector),axis=1)
	
	#tau is defined in Gagne et al. (2017, ApJS, X, Y, equation 5)
	TAU_vector = np.repeat(center_vec.reshape(1,6),num_stars,axis=0)
	if psira is not None:
		TAU_vector += PHI_vector
	
	#Take scalar products in multivariate space
	OMEGA_OMEGA = scalar_set_product_multivariate(OMEGA_vector,OMEGA_vector,precision_matrix)
	GAMMA_GAMMA = scalar_set_product_multivariate(GAMMA_vector,GAMMA_vector,precision_matrix)
	OMEGA_GAMMA = scalar_set_product_multivariate(OMEGA_vector,GAMMA_vector,precision_matrix)
	OMEGA_TAU = scalar_set_product_multivariate(OMEGA_vector,TAU_vector,precision_matrix)
	GAMMA_TAU = scalar_set_product_multivariate(GAMMA_vector,TAU_vector,precision_matrix)
	TAU_TAU = scalar_set_product_multivariate(TAU_vector,TAU_vector,precision_matrix)
	
	#If radial velocity or distance measurements are given, propagate them to the relevant scalar products
	if dist_measured is not None and dist_error is not None:
		#Find where measured distances are finite
		finite_ind = np.where(np.isfinite(dist_measured) & np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(dist_error[finite_ind],1e-3)**2
			GAMMA_GAMMA[finite_ind] += 1.0/norm
			GAMMA_TAU[finite_ind] += dist_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += dist_measured[finite_ind]**2/norm
	if rv_measured is not None and rv_error is not None:
		#Find where measured RVs are finite
		finite_ind = np.where(np.isfinite(rv_measured) & np.isfinite(rv_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(rv_error[finite_ind],1e-3)**2
			OMEGA_OMEGA[finite_ind] += 1.0/norm
			OMEGA_TAU[finite_ind] += rv_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += rv_measured[finite_ind]**2/norm
	
	#Calculate the determinant of the precision matrix unless it is given as a parameter
	if precision_matrix_determinant is None:
		precision_matrix_determinant = np.linalg.det(precision_matrix)
	if precision_matrix_determinant <= 0:
		raise ValueError('The determinant of the precision matrix bust be positive and non-zero !')
	
	#Calculate optimal distance and radial velocity
	beta = (GAMMA_GAMMA - OMEGA_GAMMA**2/OMEGA_OMEGA)/2.0
	if np.nanmin(beta) < 0:
		raise ValueError('beta has an ill-defined value !')
	gamma = OMEGA_GAMMA*OMEGA_TAU/OMEGA_OMEGA - GAMMA_TAU
	dist_optimal = (np.sqrt(gamma**2+32.0*beta) - gamma) / (4.0*beta)
	rv_optimal = (4.0 - GAMMA_GAMMA*dist_optimal**2 + GAMMA_TAU*dist_optimal)/(OMEGA_GAMMA*dist_optimal)
	
	#Create arrays that contain the measured RV and distance if available, or the optimal values otherwise
	dist_optimal_or_measured = dist_optimal
	if dist_measured is not None and dist_error is not None:
		finite_ind = np.where(np.isfinite(dist_measured) & np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			dist_optimal_or_measured[finite_ind] = dist_measured[finite_ind]
	rv_optimal_or_measured = rv_optimal
	if rv_measured is not None and rv_error is not None:
		finite_ind = np.where(np.isfinite(rv_measured) & np.isfinite(rv_error))
		if np.size(finite_ind) != 0:
			rv_optimal_or_measured[finite_ind] = rv_measured[finite_ind]
	
	#Propagate proper motion measurement errors
	EX = np.zeros(num_stars)
	EY = np.zeros(num_stars)
	EZ = np.zeros(num_stars)
	(U, V, W, EU, EV, EW) = equatorial_UVW(ra,dec,pmra,pmdec,rv_optimal_or_measured,dist_optimal_or_measured,pmra_error=pmra_error,pmdec_error=pmdec_error)
	
	#Determine by how much the diagonal of the covariance matrix must be inflated to account for the measurement errors
	covariance_matrix = np.linalg.inv(precision_matrix)
	covariance_diagonal = np.diag(covariance_matrix)
	inflation_array = np.array([EX,EY,EZ,EU,EV,EW]).transpose()
	inflation_factors = 1.0 + inflation_array**2/np.repeat(covariance_diagonal.reshape(1,6),num_stars,axis=0)
	
	#Calculate how much the determinant of the covariance matrices must be inflated
	inflation_covariance_determinant = np.exp(np.sum(np.log(inflation_factors),axis=1))
	
	#Make sure that no matrix becomes unphysical
	if np.nanmin(inflation_covariance_determinant) <= 0:
		raise ValueError('At least one covariance matrix has a negative or null determinant as a consequence of the measurement errors !')
	
	#Calculate new determinants for the precision matrices
	precision_matrix_inflated_determinant = precision_matrix_determinant/inflation_covariance_determinant
	
	#Apply this to the precision matrices
	precision_matrix_inflated = matrix_set_inflation(precision_matrix, 1.0/np.sqrt(inflation_factors))
	
	#Recalculate the scalar products with new precision matrices
	OMEGA_OMEGA = scalar_set_product_multivariate_variablemetric(OMEGA_vector,OMEGA_vector,precision_matrix_inflated)
	GAMMA_GAMMA = scalar_set_product_multivariate_variablemetric(GAMMA_vector,GAMMA_vector,precision_matrix_inflated)
	OMEGA_GAMMA = scalar_set_product_multivariate_variablemetric(OMEGA_vector,GAMMA_vector,precision_matrix_inflated)
	OMEGA_TAU = scalar_set_product_multivariate_variablemetric(OMEGA_vector,TAU_vector,precision_matrix_inflated)
	GAMMA_TAU = scalar_set_product_multivariate_variablemetric(GAMMA_vector,TAU_vector,precision_matrix_inflated)
	TAU_TAU = scalar_set_product_multivariate_variablemetric(TAU_vector,TAU_vector,precision_matrix_inflated)
	
	#If radial velocity or distance measurements are given, propagate them to the relevant scalar products
	if dist_measured is not None and dist_error is not None:
		#Find where measured distances are finite
		finite_ind = np.where(np.isfinite(dist_measured) & np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(dist_error[finite_ind],1e-3)**2
			GAMMA_GAMMA[finite_ind] += 1.0/norm
			GAMMA_TAU[finite_ind] += dist_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += dist_measured[finite_ind]**2/norm
	if rv_measured is not None and rv_error is not None:
		#Find where measured RVs are finite
		finite_ind = np.where(np.isfinite(rv_measured) & np.isfinite(rv_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(rv_error[finite_ind],1e-3)**2
			OMEGA_OMEGA[finite_ind] += 1.0/norm
			OMEGA_TAU[finite_ind] += rv_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += rv_measured[finite_ind]**2/norm
	
	#Update optimal distance and radial velocity
	beta = (GAMMA_GAMMA - OMEGA_GAMMA**2/OMEGA_OMEGA)/2.0
	if np.nanmin(beta) < 0:
		raise ValueError('beta has an ill-defined value !')
	gamma = OMEGA_GAMMA*OMEGA_TAU/OMEGA_OMEGA - GAMMA_TAU
	dist_optimal = (np.sqrt(gamma**2+32.0*beta) - gamma) / (4.0*beta)
	rv_optimal = (4.0 - GAMMA_GAMMA*dist_optimal**2 + GAMMA_TAU*dist_optimal)/(OMEGA_GAMMA*dist_optimal)
	
	#Calculate error bars on the optimal distance and radial velocity
	edist_optimal = 1.0/np.sqrt(GAMMA_GAMMA)
	erv_optimal = 1.0/np.sqrt(OMEGA_OMEGA)
	
	#Calculate final quantities for ln probability
	zeta = (TAU_TAU - OMEGA_TAU**2/OMEGA_OMEGA)/2.0
	xarg = gamma/np.sqrt(2.0*beta)
	
	lnP_coeff = -0.5*np.log(OMEGA_OMEGA) - 2.5*np.log(beta) + 0.5*np.log(precision_matrix_inflated_determinant)
	lnP_part1 = xarg**2/2.0 - zeta
	lnP_part2 = np.log(np.maximum(parabolic_cylinder_f5_mod(xarg),tiny_number))
	lnP = lnP_coeff + lnP_part1 + lnP_part2
	
	#Return ln_P if only this is required
	if lnP_only:
		return lnP
	
	#Create arrays that contain the measured RV and distance if available, or the optimal values otherwise
	dist_optimal_or_measured = dist_optimal
	edist_optimal_or_measured = edist_optimal
	if dist_measured is not None and dist_error is not None:
		finite_ind = np.where(np.isfinite(dist_measured) & np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			dist_optimal_or_measured[finite_ind] = dist_measured[finite_ind]
			edist_optimal_or_measured[finite_ind] = dist_error[finite_ind]
	rv_optimal_or_measured = rv_optimal
	erv_optimal_or_measured = erv_optimal
	if rv_measured is not None and rv_error is not None:
		finite_ind = np.where(np.isfinite(rv_measured) & np.isfinite(rv_error))
		if np.size(finite_ind) != 0:
			rv_optimal_or_measured[finite_ind] = rv_measured[finite_ind]
			erv_optimal_or_measured[finite_ind] = rv_error[finite_ind]
	
	#Calculate XYZ and UVW positions at the optimal (or measured) RV and distance
	(X, Y, Z, EX, EY, EZ) = equatorial_XYZ(ra,dec,dist_optimal_or_measured,dist_error=edist_optimal_or_measured)
	(U, V, W, EU, EV, EW) = equatorial_UVW(ra,dec,pmra,pmdec,rv_optimal_or_measured,dist_optimal_or_measured,pmra_error=pmra_error,pmdec_error=pmdec_error,rv_error=erv_optimal_or_measured,dist_error=edist_optimal_or_measured)
	XYZUVW = np.array([X,Y,Z,U,V,W]).transpose()
	EXYZUVW = np.array([EX,EY,EZ,EU,EV,EW]).transpose()
	
	#Calculate the Mahalanobis distance from the optimal position to the Gaussian model
	vec = XYZUVW - TAU_vector
	mahalanobis = np.sqrt(scalar_set_product_multivariate_variablemetric(vec,vec,precision_matrix_inflated))
	
	#Calculate the XYZ (pc) and UVW (km/s) separations from the optimal position to the center of the Gaussian model
	XYZ_sep = np.sqrt(np.sum((XYZUVW[:,0:3]-TAU_vector[:,0:3])**2,axis=1))
	UVW_sep = np.sqrt(np.sum((XYZUVW[:,3:6]-TAU_vector[:,3:6])**2,axis=1))
	
	#Calculate the 3D N-sigma distances from the optimal position to the center of the Gaussian models
	XYZ_sig = np.sqrt(scalar_set_product_multivariate_variablemetric(vec[:,0:3],vec[:,0:3],precision_matrix_inflated[:,0:3,0:3]))
	UVW_sig = np.sqrt(scalar_set_product_multivariate_variablemetric(vec[:,3:6],vec[:,3:6],precision_matrix_inflated[:,3:6,3:6]))
	
	#Store the data in a pandas dataframe
	output_structure = pd.DataFrame(np.array([lnP,dist_optimal,rv_optimal,edist_optimal,erv_optimal,X,Y,Z,U,V,W,EX,EY,EZ,EU,EV,EW,XYZ_sep,UVW_sep,XYZ_sig,UVW_sig,mahalanobis]).transpose(),columns=['LN_P','D_OPT','RV_OPT','ED_OPT','ERV_OPT','X','Y','Z','U','V','W','EX','EY','EZ','EU','EV','EW','XYZ_SEP','UVW_SEP','XYZ_SIG','UVW_SIG','MAHALANOBIS'])
	
	#Return the output table
	return output_structure
	
def parabolic_cylinder_f5_mod(x):
	"""
	Calculates the real part of the "modified" Parabolic Cylinder Function D of index v=-5.
	
	The regular function D(-5,x) is equivalent to the real part of:
		from scipy.special import pbdv
		return pbdv(-5,x)
	
	And is equivalent to the mathematical expression:
		exp(x^2/4)/24 * (sqrt(pi/2)*(x^4+6*x^2+3)*erfc(x/sqrt(2)) - exp(-x^2/2)*(x^3+5*x))
		
	The modified parabolic cylinder does away with the exp(x^2/4) term to improve numerical stability, and instead returns:
		(sqrt(pi/2)*(x^4+6*x^2+3)*erfc(x/sqrt(2)) - exp(-x^2/2)*(x^3+5*x))/24
	
	"""
	
	#Define shortcuts for efficiency
	sqrt2 = np.sqrt(2.)
	sqrt_halfpi = np.sqrt(np.pi)/sqrt2
	x_over_sqrt2 = x / sqrt2
	erfc_x_over_sqrt2 = erfc(x_over_sqrt2)
	epsilon = np.exp(-x**2/2.0)
	
	#Calculate the output
	y = 1/24.0*(sqrt_halfpi*(x**4+6.*x**2+3.)*erfc_x_over_sqrt2 - epsilon*(x**3+5.*x))
	
	return y
	
def equatorial_galactic(ra,dec):
	"""Transforms equatorial coordinates (ra,dec) to Galactic coordinates (gl,gb). All inputs must be numpy arrays of the same dimension
	
		param ra: Right ascension (degrees)
		param dec: Declination (degrees)
		output (gl,gb): Tuple containing Galactic longitude and latitude (degrees)
	"""
	
	#Check for parameter consistency
	num_stars = np.size(ra)
	if np.size(dec) != num_stars:
		raise ValueError('The dimensions ra and dec do not agree. They must all be numpy arrays of the same length.')
	
	#Compute intermediate quantities
	ra_m_ra_pol = ra - ra_pol
	sin_ra = np.sin(np.radians(ra_m_ra_pol))
	cos_ra = np.cos(np.radians(ra_m_ra_pol))
	sin_dec = np.sin(np.radians(dec))
	cos_dec = np.cos(np.radians(dec))
	
	#Compute Galactic latitude
	gamma = sin_dec_pol*sin_dec + cos_dec_pol*cos_dec*cos_ra
	gb = np.degrees(np.arcsin(gamma))
	
	#Compute Galactic longitude
	x1 = cos_dec * sin_ra
	x2 = (sin_dec - sin_dec_pol*gamma)/cos_dec_pol
	gl = l_north - np.degrees(np.arctan2(x1,x2))
	gl = (gl+360.)%(360.)
	#gl = np.mod(gl,360.0) might be better
	
	#Return Galactic coordinates tuple
	return (gl, gb)
	
def matrix_set_product_A_single(A,B):
	"""Performs matrix multiplication A#B where B is a set of N matrices. This function is more performant than looping over the N matrices if N is much larger than the matrix dimension D. A and the individual Bs must be square. The columns of A are multiplied by the rows of Bs. In IDL this function is called matrix_multiply_square_act.
	
	param A: DxD numpy array
	param B: NxDxD numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = A.shape[0]
	set_size = B.shape[0]
	if A.shape[1] != matrix_dim or B.shape[1] != matrix_dim or B.shape[2] != matrix_dim:
		raise ValueError('The dimensions D of matrices A and B do not agree - A must have dimension DxD and B must have dimension NxDxD')
	if np.size(A.shape) != 2 or np.size(B.shape) != 3:
		raise ValueError('The number of dimensions of matrices A and B are not valid - A must have dimension DxD and B must have dimension NxDxD')
	
	#Initiate resulting matrix C and perform by-element matrix multiplication
	C = np.zeros([set_size,matrix_dim,matrix_dim])
	for i in range(0,matrix_dim):
		for j in range(0,matrix_dim):
			for k in range(0,matrix_dim):
				C[:,i,j] += A[i,k] * B[:,k,j]
	
	#Return the resulting matrix
	return C
	
def matrix_vector_set_product_v_single(A,v):
	"""Performs matrix-vector multiplication A#v where A is a set of matrices and v is a single vector. This function is more performant than looping over the N sets if N is much larger than the matrix-vector dimension D. A must be square. Each column of A is multiplied by the vector v in a scalar product. In IDL this function is called matrix_vector_product_vct.
	
	param A: NxDxD numpy array
	param v: D numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = A.shape[1]
	set_size = A.shape[0]
	if A.shape[2] != matrix_dim or v.shape[0] != matrix_dim:
		raise ValueError('The dimensions D of matrix A and vector v do not agree - A must have dimension NxDxD and v must have dimension D')
	if np.size(A.shape) != 3 or np.size(v.shape) != 1:
		raise ValueError('The number of dimensions of matrix A vector v are not valid - A must have dimension NxDxD and v must have dimension D')
	
	#Initiate resulting vector w and perform by-element matrix-vector multiplication
	w = np.zeros([set_size,matrix_dim])
	for i in range(0,matrix_dim):
		for k in range(0,matrix_dim):
			w[:,i] += A[:,i,k] * v[k]
	
	#Return the resulting vector
	return w

def matrix_vector_set_product(A,v):
	"""
	Performs matrix-vector multiplication A#v where both A and v are sets of N matrices and N vectors. This function is more performant than looping over the N sets if N is much larger than the matrix-vector dimension D. A must be square. Each column of A is multiplied by the vector v in a scalar product. In IDL this function is called matrix_vector_product.
	
	param A: NxDxD numpy array
	param v: NxD numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = A.shape[1]
	set_size = A.shape[0]
	if A.shape[2] != matrix_dim or v.shape[1] != matrix_dim:
		raise ValueError('The dimensions D of matrix A and vector v do not agree - A must have dimension NxDxD and v must have dimension NxD')
	if np.size(A.shape) != 3 or np.size(v.shape) != 2:
		raise ValueError('The number of dimensions of matrix A vector v are not valid - A must have dimension NxDxD and v must have dimension NxD')
	
	#Initiate resulting vector w and perform by-element matrix-vector multiplication
	w = np.zeros([set_size,matrix_dim])
	for i in range(0,matrix_dim):
		for k in range(0,matrix_dim):
			w[:,i] += A[:,i,k] * v[:,k]
	
	#Return the resulting vector
	return w
	
def scalar_set_product_multivariate(u,v,metric):
	"""
	Performs scalar multiplication in a non-Euclidian metric u#(metric)#v. Both u and v are sets of N vectors. This function is more performant than looping over the N vectors if N is much larger than the vector dimension D. In IDL this function is called inner_product_multi.
	
	param u: NxD numpy array
	param v: NxD numpy array
	param metric: DxD numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = u.shape[1]
	set_size = u.shape[0]
	if v.shape[0] != set_size or v.shape[1] != matrix_dim:
		raise ValueError('The dimensions of vectors u and v do not agree - both must have dimension NxD')
	if metric.shape[0] != matrix_dim or metric.shape[1] != matrix_dim:
		raise ValueError('The dimensions of the metric are incompatible with vectors u and v - It must have dimension DxD where u and v have dimensions NxD')
	if np.size(u.shape) != 2 or np.size(v.shape) != 2 or np.size(metric.shape) != 2:
		raise ValueError('The number of dimensions of vectors u, v and metric matrix are not valid - u and v must have dimension NxD and metric must have dimension DxD')
	
	#Initiate resulting scalar w and perform by-element matrix multiplication
	w = np.zeros(set_size)
	for i in range(0,matrix_dim):
		for j in range(0,matrix_dim):
			w += u[:,i] * v[:,j] * metric[i,j]
	
	#Return the resulting scalar
	return w
	
def scalar_set_product_multivariate_variablemetric(u,v,metric):
	"""
	Performs scalar multiplication in a non-Euclidian metric u#(metric)#v. Both u and v are sets of N vectors, and "metric" is a set of matrices. This function is more performant than looping over the N vectors if N is much larger than the vector dimension D. In IDL this function is called inner_product_multi.
	
	param u: NxD numpy array
	param v: NxD numpy array
	param metric: NxDxD numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = u.shape[1]
	set_size = u.shape[0]
	if v.shape[0] != set_size or v.shape[1] != matrix_dim:
		raise ValueError('The dimensions of vectors u and v do not agree - both must have dimension NxD')
	if metric.shape[0] != set_size or metric.shape[1] != matrix_dim or metric.shape[2] != matrix_dim:
		raise ValueError('The dimensions of the metric are incompatible with vectors u and v - It must have dimension NxDxD where u and v have dimensions NxD')
	if np.size(u.shape) != 2 or np.size(v.shape) != 2 or np.size(metric.shape) != 3:
		raise ValueError('The number of dimensions of vectors u, v and metric matrix are not valid - u and v must have dimension NxD and metric must have dimension NxDxD')
	
	#Initiate resulting scalar w and perform by-element matrix multiplication
	w = np.zeros(set_size)
	for i in range(0,matrix_dim):
		for j in range(0,matrix_dim):
			w += u[:,i] * v[:,j] * metric[:,i,j]
	
	#Return the resulting scalar
	return w
	
def matrix_set_inflation(A,v):
	"""
	Performs the inflation of the diagonal of a single matrix A with set of factors v. This is the equivalent of v#A#v.
	
	param A: DxD numpy array
	param v: NxD numpy array
	"""
	
	#Verify matrix dimensions
	matrix_dim = A.shape[0]
	set_size = v.shape[0]
	if A.shape[1] != matrix_dim or v.shape[1] != matrix_dim:
		raise ValueError('The dimensions of matrix A vector v do not agree - A must have dimension DxD and v must have dimension NxD')
	if np.size(A.shape) != 2 or np.size(v.shape) != 2:
		raise ValueError('The number of dimensions of matrix A or vector v are not valid - A must have dimension DxD and v must have dimension NxD')
	
	#Calculate B = A#v
	B = np.empty((set_size,matrix_dim,matrix_dim))
	for i in range(0,matrix_dim):
		for j in range(0,matrix_dim):
			B[:,i,j] = A[i,j] * v[:,j]
	
	#Calculate C = v#B = v#A#v
	C = np.empty((set_size,matrix_dim,matrix_dim))
	for i in range(0,matrix_dim):
		for j in range(0,matrix_dim):
			C[:,i,j] = v[:,i] * B[:,i,j]
	
	#Return the resulting set of matrices
	return C
	
def equatorial_XYZ(ra,dec,dist,dist_error=None):
	"""
	Transforms equatorial coordinates (ra,dec) and distance to Galactic position XYZ. All inputs must be numpy arrays of the same dimension.
	
	param ra: Right ascension (degrees)
	param dec: Declination (degrees)
	param dist: Distance (parsec)
	param dist_error: Error on distance (parsec)
	
	output (X,Y,Z): Tuple containing Galactic position XYZ (parsec)
	output (X,Y,Z,EX,EY,EZ): Tuple containing Galactic position XYZ and their measurement errors, used if any measurement errors are given as inputs (parsec)
	"""
	
	#Verify keywords
	num_stars = np.size(ra)
	if np.size(dec) != num_stars or np.size(dist) != num_stars:
		raise ValueError('ra, dec and distance must all be numpy arrays of the same size !')
	if dist_error is not None and np.size(dist_error) != num_stars:
		raise ValueError('dist_error must be a numpy array of the same size as ra !')
	
	#Compute Galactic coordinates
	(gl, gb) = equatorial_galactic(ra,dec)
	
	cos_gl = np.cos(np.radians(gl))
	cos_gb = np.cos(np.radians(gb))
	sin_gl = np.sin(np.radians(gl))
	sin_gb = np.sin(np.radians(gb))
	
	X = cos_gb * cos_gl * dist
	Y = cos_gb * sin_gl * dist
	Z = sin_gb * dist
	
	if dist_error is not None:
		#X_gb = sin_gb * cos_gl * dist * np.pi/180.
		#X_gl = cos_gb * sin_gl * dist * np.pi/180.
		X_dist = cos_gb * cos_gl
		EX = np.abs(X_dist * dist_error)
		Y_dist = cos_gb * sin_gl
		EY = np.abs(Y_dist * dist_error)
		Z_dist = sin_gb
		EZ = np.abs(Z_dist * dist_error)
		return (X, Y, Z, EX, EY, EZ)
	else:
		return (X, Y, Z)
		
def equatorial_UVW(ra,dec,pmra,pmdec,rv,dist,pmra_error=None,pmdec_error=None,rv_error=None,dist_error=None):
	"""
	Transforms equatorial coordinates (ra,dec), proper motion (pmra,pmdec), radial velocity and distance to space velocities UVW. All inputs must be numpy arrays of the same dimension.
	
	param ra: Right ascension (degrees)
	param dec: Declination (degrees)
	param pmra: Proper motion in right ascension (milliarcsecond per year). 	Must include the cos(delta) term
	param pmdec: Proper motion in declination (milliarcsecond per year)
	param rv: Radial velocity (kilometers per second)
	param dist: Distance (parsec)
	param ra_error: Error on right ascension (degrees)
	param dec_error: Error on declination (degrees)
	param pmra_error: Error on proper motion in right ascension (milliarcsecond per year)
	param pmdec_error: Error on proper motion in declination (milliarcsecond per year)
	param rv_error: Error on radial velocity (kilometers per second)
	param dist_error: Error on distance (parsec)
	
	output (U,V,W): Tuple containing Space velocities UVW (kilometers per second)
	output (U,V,W,EU,EV,EW): Tuple containing Space velocities UVW and their measurement errors, used if any measurement errors are given as inputs (kilometers per second)
	"""
	
	#Verify keywords
	num_stars = np.size(ra)
	if np.size(dec) != num_stars or np.size(pmra) != num_stars or np.size(pmdec) != num_stars or np.size(dist) != num_stars:
		raise ValueError('ra, dec, pmra, pmdec, rv and distance must all be numpy arrays of the same size !')
	if pmra_error is not None and np.size(pmra_error) != num_stars:
		raise ValueError('pmra_error must be a numpy array of the same size as ra !')
	if pmdec_error is not None and np.size(pmdec_error) != num_stars:
		raise ValueError('pmdec_error must be a numpy array of the same size as ra !')
	if rv_error is not None and np.size(rv_error) != num_stars:
		raise ValueError('rv_error must be a numpy array of the same size as ra !')
	if dist_error is not None and np.size(dist_error) != num_stars:
		raise ValueError('dist_error must be a numpy array of the same size as ra !')
	
	#Compute elements of the T matrix
	cos_ra = np.cos(np.radians(ra))
	cos_dec = np.cos(np.radians(dec))
	sin_ra = np.sin(np.radians(ra))
	sin_dec = np.sin(np.radians(dec))
	T1 = TGAL[0,0]*cos_ra*cos_dec + TGAL[0,1]*sin_ra*cos_dec + TGAL[0,2]*sin_dec
	T2 = -TGAL[0,0]*sin_ra + TGAL[0,1]*cos_ra
	T3 = -TGAL[0,0]*cos_ra*sin_dec - TGAL[0,1]*sin_ra*sin_dec + TGAL[0,2]*cos_dec
	T4 = TGAL[1,0]*cos_ra*cos_dec + TGAL[1,1]*sin_ra*cos_dec + TGAL[1,2]*sin_dec
	T5 = -TGAL[1,0]*sin_ra + TGAL[1,1]*cos_ra
	T6 = -TGAL[1,0]*cos_ra*sin_dec - TGAL[1,1]*sin_ra*sin_dec + TGAL[1,2]*cos_dec
	T7 = TGAL[2,0]*cos_ra*cos_dec + TGAL[2,1]*sin_ra*cos_dec + TGAL[2,2]*sin_dec
	T8 = -TGAL[2,0]*sin_ra + TGAL[2,1]*cos_ra
	T9 = -TGAL[2,0]*cos_ra*sin_dec - TGAL[2,1]*sin_ra*sin_dec + TGAL[2,2]*cos_dec
	
	#Calculate UVW
	reduced_dist = kappa*dist
	U = T1*rv + T2*pmra*reduced_dist + T3*pmdec*reduced_dist
	V = T4*rv + T5*pmra*reduced_dist + T6*pmdec*reduced_dist
	W = T7*rv + T8*pmra*reduced_dist + T9*pmdec*reduced_dist
	
	#Return only (U, V, W) tuple if no errors are set
	if pmra_error is None and pmdec_error is None and rv_error is None and dist_error is None:
		return (U, V, W)
		
	#Propagate errors if they are specified
	if pmra_error is None:
		pmra_error = np.zeros(num_stars)
	if pmdec_error is None:
		pmdec_error = np.zeros(num_stars)
	if rv_error is None:
		rv_error = np.zeros(num_stars)
	if dist_error is None:
		dist_error = np.zeros(num_stars)
	reduced_dist_error = kappa*dist_error
	
	#Calculate derivatives
	T23_pm = np.sqrt((T2*pmra)**2+(T3*pmdec)**2)
	T23_pm_error = np.sqrt((T2*pmra_error)**2+(T3*pmdec_error)**2)
	EU_rv = T1 * rv_error
	EU_pm = T23_pm_error * reduced_dist
	EU_dist = T23_pm * reduced_dist_error
	EU_dist_pm = T23_pm_error * reduced_dist_error
	
	T56_pm = np.sqrt((T5*pmra)**2+(T6*pmdec)**2)
	T56_pm_error = np.sqrt((T5*pmra_error)**2+(T6*pmdec_error)**2)
	EV_rv = T4 * rv_error
	EV_pm = T56_pm_error * reduced_dist
	EV_dist = T56_pm * reduced_dist_error
	EV_dist_pm = T56_pm_error * reduced_dist_error

	T89_pm = np.sqrt((T8*pmra)**2+(T9*pmdec)**2)
	T89_pm_error = np.sqrt((T8*pmra_error)**2+(T9*pmdec_error)**2)
	EW_rv = T7 * rv_error
	EW_pm = T89_pm_error * reduced_dist
	EW_dist = T89_pm * reduced_dist_error
	EW_dist_pm = T89_pm_error * reduced_dist_error
	
	#Calculate error bars
	EU = np.sqrt(EU_rv**2 + EU_pm**2 + EU_dist**2 + EU_dist_pm**2)
	EV = np.sqrt(EV_rv**2 + EV_pm**2 + EV_dist**2 + EV_dist_pm**2)
	EW = np.sqrt(EW_rv**2 + EW_pm**2 + EW_dist**2 + EW_dist_pm**2)
	
	#Return measurements and error bars
	return (U, V, W, EU, EV, EW)
	