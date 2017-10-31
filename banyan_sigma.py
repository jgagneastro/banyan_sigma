"""
 NAME:
       BANYAN_SIGMA

 PURPOSE:
       Calculate the membership probability that a given astrophysical object belong to one of the currently
       known 27 young associations within 150 pc of the Sun, using Bayesian inference. This tool uses the sky
       position and proper motion measurements of an object, with optional radial velocity (RV) and distance (D)
       measurements, to derive a Bayesian membership probability. By default, the priors are adjusted such that
       a probability treshold of 90% will recover 50%, 68%, 82% or 90% of true association members depending on
       what observables are input (only sky position and proper motion, with RV, with D, with both RV and D,
       respectively).
       
       Please see Gagné et al. 2017 (ApJ, XX, XX; Arxiv YY, YY) for more detail.
       
       An online version of this tool is available for 1-object queries at http://www.astro.umontreal.ca/~gagne/banyansigma.php.
       
 REQUIREMENTS:
       (1) A fits file containing the parameters of the multivariate Gaussian models of each Bayesian hypothesis must be included
           at /data/banyan_sigma_parameters.fits in the directory where BANYAN_SIGMA() is compiled. 
           The fits file can be written with the IDL MWRFITS.PRO function from an IDL array of structures of N elements, where N is the total
           number of multivariate Gaussians used in the models of all Bayesian hypotheses. Each element of this structure contains
           the following information:
           NAME - The name of the model (scalar string).
           CENTER_VEC - Central XYZUVW position of the model (6D vector, in units of pc and km/s).
           COVARIANCE_MATRIX - Covariance matrix in XYZUVW associated with the model (6x6 matrix, in mixed units of pc and km/s).
           PRECISION_MATRIX - (Optional) Matrix inverse of COVARIANCE_MATRIX, to avoid re-calculating it many times (6x6 matrix).
           LN_NOBJ - (Optional) Natural logarithm of the number of objects used to build the synthetic model (scalar). This is not used in banyan_sigma().
           COVARIANCE_DETERM - (Optional) Determinant of the covariance matrix, to avoid re-calculating it many times (scalar).
           PRECISION_DETERM - (Optional) Determinant of the precision matrix, to avoid re-calculating it many times (scalar).
           LN_ALPHA_K - (Optional) Natural logarithm of the alpha_k inflation factors that ensured a fixed rate of true positives
                        at a given Bayesian probability treshold. See Gagné et al. 2017 (ApJ, XX, XX) for more detail (scalar or 4-elements vector). This is not used in BANYAN_SIGMA.
           LN_PRIOR - Natural logarithm of the Bayesian prior (scalar of 4-elements vector). When this is a 4-elements vector,
                      the cases with only proper motion, proper motion + radial velocity, proper motion + distance or proper motion + radial velocity + distance
                      will be used with the corresponding element of the LN_PRIOR vector.
           LN_PRIOR_OBSERVABLES - Scalar string or 4-elements vector describing the observing modes used for each element of ln_prior.
                                  This is not used in banyan_sigma().
           COEFFICIENT - Coefficient (or weight) for multivariate Gaussian mixture models. This will only be used if more than one element of the
                         parameters array have the same model name (see below).  
           
           In Python, this fits file is read with the Astropy.Tables routine.
           
           When more than one elements have the same model name, BANYAN_SIGMA will use the COEFFICIENTs to merge its Bayesian probability,
           therefore representing the hypothesis with a multivariate Gaussian model mixture.
           
       (2) (Optional) A fits file containing the various performance metrics (true positive rate, false positive rate, positive
           predictive value) as a function of the Bayesian probability treshold, for each young association. Each element of this structure contains
           the following information:
           NAME - The name of the model (scalar string).
           PROBS - N-elements array containing a list of Bayesian probabilities (%).
           TPR - Nx4-elements array containing the rate of true positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
           FPR - Nx4-elements array containing the rate  of false positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
           MCC - Nx4-elements array containing the Matthews Correlation Coefficients that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
           PPV - Nx4-elements array containing the Positive Predictive Values that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
           
           Each component of the 4-elements dimension of TPR, FPR, MCC and PPV corresponds to a different mode of input data,
           see the description of "LN_PRIOR" above for more detail.
           
           When this fits file is used, the Bayesian probabilities of each star will be associated with a TPR, FPR, MCC and PPV values in the METRICS sub-structure of
           the output structure.
           
           This file must be located at /data/banyan_sigma_metrics.fits in the directory where BANYAN_SIGMA.pro is compiled.
           
 CALLING SEQUENCE:
 
       OUTPUT_STRUCTURE = BANYAN_SIGMA(stars_data=None, column_names=None, hypotheses=None, ln_priors=None, ntargets_max=None, 
           ra=None, dec=None, pmra=None, pmdec=None, epmra=None, epmdec=None, dist=None, edist=None, rv=None, erv=None,
           psira=None, psidec=None, epsira=None, epsidec=None, plx=None, eplx=None,
           constraint_dist_per_hyp=None, constraint_edist_per_hyp=None,
          unit_priors=True/False, lnp_only=True/False, no_xyz=True/False, use_rv=True/False, use_dist=True/False, 
          use_plx=True/False, use_psi=True/False)

 OPTIONAL INPUTS:
       stars_data - An structured array that contains at least the following informations:
                    ra, dec, pmra, pmdec, epmra, and epmdec. It can also optionally contain the informations on
                    rv, erv, dist, edist, plx, eplx, psira, psidec, epsira, epsidec. See the corresponding keyword
                    descriptions for more information.
                    If this input is not used, the keywords ra, dec, pmra, pmdec, epmra, and epmdec must all be specified.
       column_names - A Python dictionary that contains the names of the "stars_data" columns columns which differ from the
                     default values listed above. For example, column_names = {'RA':'ICRS_RA'} can be used to specify that
                     the RA values are listed in the column of stars_data named ICRS_RA.
       ra - Right ascension (decimal degrees). A N-elements array can be specified to calculate the Bayesian probability of several
            stars at once, but then all mandatory inputs must also be N-elements arrays.
       dec - Declination (decimal degrees).
       pmra - Proper motion in the right ascension direction (mas/yr, must include the cos(dec) factor).
       pmdec - Proper motion in the declination direction (mas/yr).
       epmra - Measurement error on the proper motion in the right ascension direction (mas/yr, must not include the cos(dec) factor).
       epmdec -  Measurement error on the proper motion in the declination direction (mas/yr).
       rv - Radial velocity measurement to be included in the Bayesian probability (km/s).
            If this keyword is set, erv must also be set.
            A N-elements array must be used if N stars are analyzed at once.
       erv - Measurement error on the radial velocity to be included in the Bayesian probability (km/s).
             A N-elements array must be used if N stars are analyzed at once.
       dist - Distance measurement to be included in the Bayesian probability (pc).
              By default, the banyan_sigma() Bayesian priors are meant for this keyword to be used with trigonometric distances only.
              Otherwise, the rate of true positives may be far from the nominal values described in Gagné et al. (ApJS, 2017, XX, XX).
              If this keyword is set, edist must also be set.
              A N-elements array must be used if N stars are analyzed at once.
       edist - Measurement error on the distance to be included in the Bayesian probability (pc).
               A N-elements array must be used if N stars are analyzed at once.
       plx - Parallax measurement to be included in the Bayesian probability (mas). The distance will be approximated with dist = 1000/plx.
             If this keyword is set, eplx must also be set.
             A N-elements array must be used if N stars are analyzed at once.
       eplx - Measurement error on the parallax to be included in the Bayesian probability (mas).
              The distance error will be approximated with edist = 1000/plx**2*eplx.
              A N-elements array must be used if N stars are analyzed at once.
       psira - Parallax motion factor PSIRA described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr.
                If this keyword is set, the corresponding psidec, epsira and epsidec keywords must also be set.
                This measurement is only useful when proper motions are estimated from two single-epoch astrometric
                measurements. It captures the dependence of parallax motion as a function of distance, and allows
                banyan_sigma() to shift the UVW center of the moving group models, which is equivalent to
                correctly treating the input "proper motion" pmra, pmdec, epmra, epmdec as a true apparent motion.
                This keyword should *not* be used if proper motions were derived from more than two epochs, or if
                they were obtained from a full parallax solution.
                A N-elements array must be used if N stars are analyzed at once.
       psidec - Parallax motion factor psidec described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr.
                 A N-elements array must be used if N stars are analyzed at once.
       epsira - Measurement error on the parallax motion factor psira described in Gagné et al. (ApJS, 2017, XX, XX),
                 in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
       epsidec - Measurement error on the parallax motion factor psidec described in Gagné et al. (ApJS, 2017, XX, XX),
                  in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
       ntargets_max - (default 10^6). Maximum number of objects to run at once in BANYAN_SIGMA to avoid saturating the RAM.
                      If more targets are supplied, banyan_sigma() is run over a loop of several batches of ntargets_max objects. 
       hypotheses - The list of Bayesian hypotheses to be considered. They must all be present in the parameters fits file
                    (See REQUIREMENTS #1 above).
       ln_priors - An dictionary that contains the natural logarithm of Bayesian priors that should be *multiplied with the
                   default priors* (use unit_priors=True if you want only ln_priors to be considered). The structure must contain the name
                   of each hypothesis as keys, and the associated scalar value of the natural logarithm of the Bayesian prior for each key. 
       constraint_dist_per_hyp - A structured array that contains a distance constraint (in pc).
                   Each of the Bayesian hypotheses must be included as keys and the distance must be specified as its
                   associated scalar value. constraint_edist_per_hyp must also be specified if constraint_dist_per_hyp is specified.
                   This keyword is useful for including spectro-photometric distance constraints that depend on the age of the young association or field.
       constraint_edist_per_hyp - A structured array that contains a measurement
                   error on the distance constraint (in pc). Each of the Bayesian hypotheses must be included as keys and the
                   distance error must be specified as its associated scalar value.  

 OPTIONAL INPUT KEYWORD:
       unit_priors - If this keyword is set, all default priors are set to 1 (but they are still overrided by manual priors input with the keyword ln_priors).
       lnp_only - If this keyword is set, only Bayesian probabilities will be calculated and returned.
       no_xyz - If this keyword is set, the width of the spatial components of the multivariate Gaussian will be widened by a large
                 factor, so that the XYZ components are effectively ignored. This keyword must be used with extreme caution as it will
                 generate a significant number of false-positives and confusion between the young associations.
       use_rv - Use any radial velocity values found in the stars_data input structure.
       use_dist - Use any distance values found in the stars_data input structure.
       use_plx - Use any parallax values found in the stars_data input structure.
       use_psi - Use any psira, psidec values found in the stars_data input structure.

 OUTPUT:
      This routine outputs a structured array, with the following keys:
      NAME - The name of the object (as taken from the input structure).
      ALL - A structure that contains the Bayesian probability (0 to 1) for each of the associations (as individual keys).
      METRICS - A structure that contains the performance metrics associated with the global Bayesian probability of this target.
                This sub-structure contains the following keys:
        TPR - Rate of true positives expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
        FPR - Rate of false positives (from the field) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
        PPV - Positive Predictive Value (sample contamination) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
      [ASSOCIATION_1] - Sub-structure containing the relevant details for assiciation [ASSOCIATION_1].
      [ASSOCIATION_2] - (...).
      (...).
      [ASSOCIATION_N] - (...).
                        These sub-structures contain the following keys:
        HYPOTHESIS - Name of the association
        PROB - Bayesian probability (0 to 1)
        D_OPT - Optimal distance (pc) that maximizes the Bayesian likelihood for this hypothesis.
        RV_OPT - Optimal radial velocity (km/s) that maximizes the Bayesian likelihood for this hypothesis.
        ED_OPT - Error on the optimal distance (pc), which approximates the 68% width of how the likelihood varies with distance.
        ERV_OPT - Error on the optimal radial velocity (km/s), which approximates the 68% width of how the likelihood varies with
                  radial velocity.
        XYZUVW - 6-dimensional array containing the XYZ and UVW position of the star at the measured radial velocity and/or
                 distance, or the optimal radial velocity and/or distance when the first are not available (units of pc and km/s).
        EXYZUVW - Errors on XYZUVW (units of pc and km/s).
        XYZ_SEP - Separation between the optimal or measured XYZ position of the star and the center of the multivariate Gaussian
                  model of this Bayesian hypothesis (pc).
        UVW_SEP - Separation between the optimal or measured UVW position of the star and the center of the multivariate Gaussian
                  model of this Bayesian hypothesis (km/s).
        XYZ_SEP - N-sigma separation between the optimal or measured XYZ position of the star and the multivariate Gaussian model
                  of this Bayesian hypothesis (no units).
        UVW_SEP - N-sigma separation between the optimal or measured UVW position of the star and the multivariate Gaussian model
                  of this Bayesian hypothesis (no units).
        MAHALANOBIS - Mahalanobis distance between the optimal or measured XYZUVW position of the star and the multivariate Gaussian
                      model. A Mahalanobis distance is a generalization of a 6D N-sigma distance that accounts for covariances. 
      BESTYA_STR - A sub-structure similar to those described above for the most probable young association (ignoring the field possibility).
      YA_PROB - The Bayesian probability (0 to 1) that this object belongs to any young association (i.e., excluding the field).
      LIST_PROB_YAS - A list of young associations with at least 5% Bayesian probability. Their relative probabilities (%) are specified
                      between parentheses.
      BEST_HYP - Most probable Bayesian hypothesis (including the field)
      BEST_YA - Most probable single young association.

 MODIFICATION HISTORY:
       WRITTEN, Olivier Loubier, July, 12 2017
       MODIFIED, Jonathan Gagne, October, 25 2017
         Added several options, comments and header, performance and syntax modifications.
"""

#Import the necessary packages
import numpy as np #Numpy maths
from scipy.special import erfc
import os #Access to environment variables
from astropy.table import Table #Reading astro-formatted tables
import warnings #Raise user-defined Python warnings
import pdb #Debugging
from scipy.stats import describe #Useful for debugging

#Initiate some global constants
kappa = 0.004743717361 #1 AU/yr to km/s divided by 1000.
#For some reason "from astropy import units as u; kappa=u.au.to(u.km)/u.year.to(u.s)" is far less precise

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
def banyan_sigma(stars_data=None,column_names=None,hypotheses=None,ln_priors=None,ntargets_max=1e6,ra=None,dec=None,pmra=None,pmdec=None,epmra=None,epmdec=None,dist=None,edist=None,rv=None,erv=None,psira=None,psidec=None,epsira=None,epsidec=None,plx=None,eplx=None,constraint_dist_per_hyp=None,constraint_edist_per_hyp=None,unit_priors=False,lnp_only=False,no_xyz=False,use_rv=False,use_dist=False,use_plx=False,use_psi=False):
	
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
	
	#Check if a column named PLX, DIST, RV, PSIRA, etc. exist in stars_data but not in column_names. If this is the case, issue a warning so that the user understands that some data are not being considered.
	if stars_data is not None:
		if 'PLX' in stars_data.keys() and 'PLX' not in column_names.keys():
			warnings.warn('Parallaxes (PLX) were not read from the input data, because the PLX key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with the use_plx=True keyword to read them.')
		if 'DIST' in stars_data.keys() and 'DIST' not in column_names.keys():
			warnings.warn('Distances (DIST) were not read from the input data, because the DIST key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with the use_dist=True keyword to read them.')
		if 'RV' in stars_data.keys() and 'RV' not in column_names.keys():
			warnings.warn('Radial velocities (RV) were not read from the input data, because the RV key was not included in the column_names keyword of banyan_sigma(). You can also call banyan_sigma() with use_rv=True to read them.')
		if ('PSIRA' in stars_data.keys() and 'PSIRA' not in column_names.keys()) or ('PSIDEC' in stars_data.keys() and 'PSIDEC' not in column_names.keys()):
			warnings.warn('The PSI parameters (PSIRA,PSIDEC) were not read from the input data, because the PSIRA and PSIDEC keys were not included in the column_data keyword of banyan_sigma(). You can also call banyan_sigma() with use_psi=True keyword to read them.')
	
	#Create a table of data for BANYAN SIGMA to use
	if ra is not None:
		nobj = np.size(ra)
		zeros = np.zeros(nobj)
		rows = np.array([ra,dec,pmra,pmdec,epmra,epmdec,zeros,zeros,zeros,zeros])
	if ra is None:
		nobj = np.size(stars_data)
		zeros = np.zeros(nobj)
		rows = np.array([stars_data[column_names['RA']],stars_data[column_names['DEC']],stars_data[column_names['PMRA']],stars_data[column_names['PMDEC']],stars_data[column_names['EPMRA']],stars_data[column_names['EPMDEC']],zeros,zeros,zeros,zeros])
	data_table = Table(rows=rows.transpose(),names=('RA','DEC','PMRA','PMDEC','EPMRA','EPMDEC','PSIRA','PSIDEC','EPSIRA','EPSIDEC'))
	
	#Fill up the data table with stars_data if it is specified
	for keys in column_names.keys():
		if (keys == 'PLX') or (keys == 'EPLX'):
			continue
		data_table[keys] = stars_data[column_names[keys]]
	if 'PLX' in column_names.keys():
		data_table['DIST'] = 1e3/stars_data[column_names['PLX']]
	if 'PLX' in column_names.keys() and 'EPLX' in column_names.keys():
		data_table['EDIST'] = 1e3/stars_data[column_names['PLX']]**2*stars_data[column_names['EPLX']]
	
	#Transform parallaxes to distances directly in data_table
	if 'PLX' in data_table.keys() and 'EPLX' in data_table.keys():
		data_table['EDIST'] = 1e3/data_table['PLX']**2*data_table['EPLX']
		data_table.remove_column('EPLX')
	if 'PLX' in data_table.keys():
		data_table['DIST'] = 1e3/data_table['PLX']
		data_table.remove_column('PLX')
	
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
	
	#Data file containing the parameters of Bayesian hypotheses
	parameters_file = os.path.dirname(__file__)+os.sep+'data'+os.sep+'banyan_sigma_parameters.fits'
	
	#Check if the file exists
	if not os.path.isfile(parameters_file):
		raise ValueError('The multivariate Gaussian parameters file could not be found ! Please make sure that you did not move "'+path_sep()+'data'+path_sep()+'banyan_sigma_parameters.fits" from the same path as the Python file banyan_sigma.py !')
	
	#Read the parameters of Bayesian hypotheses
	parameters_str = Table.read(parameters_file,format='fits')
	#Remove white spaces in names
	parameters_str['NAME'] = np.chararray.strip(np.array(parameters_str['NAME']))
	npar = np.size(parameters_str)
	
	#Fix names that start with a number
	for i in range(parameters_str['NAME']):
		if parameters_str['NAME'][i][0].isdigit():
			parameters_str['NAME'][i] = '_'+parameters_str['NAME'][i]
			
	[word[0].isdigit() for word in parameters_str['NAME']]
	
	pdb.set_trace()
	
	#Issue error if measurement errors are not listed
	if 'PLX' in column_names.keys() and 'EPLX' not in column_names.keys():
		raise ValueError('If parallaxes (PLX) are specified, error bars (EPLX) must also be specified')
		
	if 'DIST' in column_names.keys() and 'EDIST' not in column_names.keys():
		raise ValueError('If distances (DIST) are specified, error bars (EDIST) must also be specified')
	
	if 'RV' in column_names.keys() and 'ERV' not in column_names.keys():
		raise ValueError('If radial velocities (RV) are specified, error bars (ERV) must also be specified')
	
	#Issue an error if both DIST and PLX are given in column_names
	if 'PLX' in column_names.keys() and 'DIST' in column_names.keys():
		raise ValueError('Distances (DIST) and parallaxes (PLX) cannot both be specified')
	
	#Issue an error if only parts of the PSI measurements or errors are given
	if 'PSIRA' in column_names.keys() ^ 'PSIDEC' in column_names.keys():
		raise ValueError('If the PSIRA keyword is specified, PSIDEC must also be specified')
	if ('PSIRA' in column_names.keys() and 'EPSIRA' not in column_names.keys()) or ('PSIDEC' in column_names.keys() and 'EPSIDEC' not in column_names.keys()):
		raise ValueError('If the PSIRA and PSIDEC keywords are specified, EPSIRA and EPSIDEC must also be specified')
	
	#Issue an error if some of the column_names entries are not present in stars_data
	for key in column_names.keys():
		if column_names[key] not in stars_data.keys():
			raise ValueError('The "'+key+'" keyword is not present in the input data structure')
	
	#Read the multivariate Gaussian parameters 
	params_file_path = os.environ['BANYAN_SIGMA_PARAMETERS']
	if not os.path.isfile(params_file_path):
		raise ValueError('The Gaussian parameters file could not be found at "'+params_file_path+'" !')
	parameters_str = Table.read(params_file_path,format='fits')
	
	#Read distance measurements from either the PLX or DIST keywords
	dist_measured = None
	dist_error = None
	if 'PLX' in column_names.keys():
		plx_measured = stars_data[column_names['PLX']]
		plx_error = stars_data[column_names['EPLX']]
		dist_measured = 1e3/plx_measured
		dist_error = 1e3/plx_measured**2*plx_error
	if 'DIST' in column_names.keys():
		dist_measured = stars_data[column_names['DIST']]
		dist_error = stars_data[column_names['EDIST']]
	
	#Read RV measurements
	rv_measured = None
	rv_error = None
	if 'RV' in column_names.keys():
		rv_measured = stars_data[column_names['RV']]
		rv_error = stars_data[column_names['ERV']]
	
	#Read PSI measurements
	psira = None
	psidec = None
	psira_error = None
	psidec_error = None
	if 'PSIRA' in column_names.keys():
		psira = stars_data[column_names['PSIRA']]
		psidec = stars_data[column_names['PSIDEC']]
		psira_error = stars_data[column_names['EPSIRA']]
		psidec_error = stars_data[column_names['EPSIDEC']]
	
	i = 0
	outi = banyan_sigma_solve_multivar(stars_data[column_names['RA']],stars_data[column_names['DEC']],stars_data[column_names['PMRA']],stars_data[column_names['PMDEC']],stars_data[column_names['EPMRA']],stars_data[column_names['EPMDEC']],rv_measured=rv_measured,rv_error=rv_error,dist_measured=dist_measured,dist_error=dist_error,psira=psira,psidec=psidec,psira_error=psira_error,psidec_error=psidec_error,precision_matrix=parameters_str[i]['PRECISION_MATRIX'],center_vec=parameters_str[i]['CENTER_VEC'],precision_matrix_determinant=parameters_str[i]['PRECISION_DETERM'])
	pdb.set_trace()
	
	lnP = 0.
	return lnP
	
def banyan_sigma_solve_multivar(ra,dec,pmra,pmdec,pmra_error,pmdec_error,precision_matrix=None,center_vec=None,rv_measured=None,rv_error=None,dist_measured=None,dist_error=None,psira=None,psidec=None,psira_error=None,psidec_error=None,lnP_only=False,precision_matrix_determinant=None):
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
	A_matrix = np.empty((num_stars,3,3))
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
		varphi_vector_sub = np.array([np.zeros(num_stars),np.array(kappa*psira), np.array(kappa*psidec)])
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
		tau_vector += PHI_vector
	
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
		finite_ind = np.where(np.isfinite(dist_measured) and np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(dist_error[finite_ind],1e-3)**2
			GAMMA_GAMMA[finite_ind] += 1.0/norm
			GAMMA_TAU[finite_ind] += dist_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += dist_measured[finite_ind]**2/norm
	if rv_measured is not None and rv_error is not None:
		#Find where measured RVs are finite
		finite_ind = np.where(np.isfinite(rv_measured) and np.isfinite(rv_error))
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
		finite_ind = np.where(np.isfinite(dist_measured) and np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			dist_optimal_or_measured[finite_ind] = dist_measured[finite_ind]
	rv_optimal_or_measured = rv_optimal
	if rv_measured is not None and rv_error is not None:
		finite_ind = np.where(np.isfinite(rv_measured) and np.isfinite(rv_error))
		if np.size(finite_ind) != 0:
			rv_optimal_or_measured[finite_ind] = rv_measured[finite_ind]
	
	#Propagate proper motion measurement errors
	#equatorial_XYZ(ra,dec,dist,dist_error=dist_error)
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
		finite_ind = np.where(np.isfinite(dist_measured) and np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			norm = np.maximum(dist_error[finite_ind],1e-3)**2
			GAMMA_GAMMA[finite_ind] += 1.0/norm
			GAMMA_TAU[finite_ind] += dist_measured[finite_ind]/norm
			TAU_TAU[finite_ind] += dist_measured[finite_ind]**2/norm
	if rv_measured is not None and rv_error is not None:
		#Find where measured RVs are finite
		finite_ind = np.where(np.isfinite(rv_measured) and np.isfinite(rv_error))
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
	lnP_part2 = np.log(np.maximum(parabolic_cylinder_f5_mod(xarg),0))
	lnP = lnP_coeff + lnP_part1 + lnP_part2
	
	#Return ln_P if only this is required
	if lnP_only:
		return lnP
	
	#Create arrays that contain the measured RV and distance if available, or the optimal values otherwise
	dist_optimal_or_measured = dist_optimal
	edist_optimal_or_measured = edist_optimal
	if dist_measured is not None and dist_error is not None:
		finite_ind = np.where(np.isfinite(dist_measured) and np.isfinite(dist_error))
		if np.size(finite_ind) != 0:
			dist_optimal_or_measured[finite_ind] = dist_measured[finite_ind]
			edist_optimal_or_measured[finite_ind] = dist_error[finite_ind]
	rv_optimal_or_measured = rv_optimal
	erv_optimal_or_measured = erv_optimal
	if rv_measured is not None and rv_error is not None:
		finite_ind = np.where(np.isfinite(rv_measured) and np.isfinite(rv_error))
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
	
	#Store relevant data in an output table
	output_structure = Table((lnP,dist_optimal,rv_optimal,edist_optimal,erv_optimal,XYZUVW,EXYZUVW,XYZ_sep,UVW_sep,XYZ_sig,UVW_sig,mahalanobis),names=('LN_P','D_OPT','RV_OPT','ED_OPT','ERV_OPT','XYZUVW','EWXYZUVW','XYZ_SEP','UVW_SEP','XYZ_SIG','UVW_SIG','MAHALANOBIS'))
	
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
	
	
	#ra_pol,dec_pol,l_north,sin_dec_pol,cos_dec_pol
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
	
	#Return galactic coordinates tuple
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
	