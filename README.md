# Banyan Σ (Python 3)
A Bayesian classifier to identify members of the 27 nearest young associations within 150 pc of the Sun.

This is the Python version of BANYAN Σ. The IDL version can be found at https://github.com/jgagneastro/banyan_sigma_idl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1165085.svg)](https://doi.org/10.5281/zenodo.1165085)

## PURPOSE:
Calculate the membership probability that a given astrophysical object belongs to one of the currently known 27 young associations within 150 pc of the Sun, using Bayesian inference. This tool uses the sky position and proper motion measurements of an object, with optional radial velocity (RV) and distance (D) measurements, to derive a Bayesian membership probability. By default, the priors are adjusted such that a probability treshold of 90% will recover 50%, 68%, 82% or 90% of true association members depending on what observables are input (only sky position and proper motion, with RV, with D, with both RV and D, respectively).
       
Please see Gagné et al. 2018 (accepted for publication in ApJS, http://adsabs.harvard.edu/abs/2018arXiv180109051G) for more detail.
       
An online version of this tool is available for 1-object queries at http://www.exoplanetes.umontreal.ca/banyan/banyansigma.php.
       
## REQUIREMENTS:
(1) This code requires Python 3 to run properly.

(2) A fits file containing the parameters of the multivariate Gaussian models of each Bayesian hypothesis must be included at /data/banyan_sigma_parameters.fits in the directory where BANYAN_SIGMA() is compiled. The file provided with this release corresponds to the set of 27 young associations described in Gagné et al. (2018). The fits file can be written with the IDL MWRFITS.PRO function from an IDL array of structures of N elements, where N is the total number of multivariate Gaussians used in the models of all Bayesian hypotheses. Each element of this structure contains the following information:

       NAME: The name of the model (scalar string).
       CENTER_VEC: Central XYZUVW position of the model (6D vector, in units of pc and km/s).
       COVARIANCE_MATRIX: Covariance matrix in XYZUVW associated with the model (6x6 matrix, in mixed units of pc and km/s).
       PRECISION_MATRIX: (Optional) Matrix inverse of COVARIANCE_MATRIX, to avoid re-calculating it many times (6x6 matrix).
       LN_NOBJ: (Optional) Natural logarithm of the number of objects used to build the synthetic model (scalar). This is not used in banyan_sigma().
       COVARIANCE_DETERM: (Optional) Determinant of the covariance matrix, to avoid re-calculating it many times (scalar).
       PRECISION_DETERM: (Optional) Determinant of the precision matrix, to avoid re-calculating it many times (scalar).
       LN_ALPHA_K: (Optional) Natural logarithm of the alpha_k inflation factors that ensured a fixed rate of true positives at a given Bayesian probability treshold. See Gagné et al. (2018) for more detail (scalar or 4-elements vector). This is not used in BANYAN_SIGMA.
       LN_PRIOR: Natural logarithm of the Bayesian prior (scalar of 4-elements vector). When this is a 4-elements vector, the cases with only proper motion, proper motion + radial velocity, proper motion + distance or proper motion + radial velocity + distance will be used with the corresponding element of the LN_PRIOR vector.
       LN_PRIOR_OBSERVABLES: Scalar string or 4-elements vector describing the observing modes used for each element of ln_prior. This is not used in banyan_sigma().
       COEFFICIENT: Coefficient (or weight) for multivariate Gaussian mixture models. This will only be used if more than one element of the parameters array have the same model name (see below).  
           
In Python, this fits file is read with the Astropy.Tables routine because it requires multi-dimensional columns. When more than one elements have the same model name, BANYAN_SIGMA will use the COEFFICIENTs to merge its Bayesian probability, therefore representing the hypothesis with a multivariate Gaussian model mixture. This is how the Galactic field is represented in Gagné et al. (2018).

(3) (Optional) A fits file containing the various performance metrics (true positive rate, false positive rate, positive predictive value) as a function of the Bayesian probability treshold, for each young association. Each element of this structure contains the following information:

       NAME: The name of the model (scalar string).
       PROBS: N-elements array containing a list of Bayesian probabilities (%).
       TPR: Nx4-elements array containing the rate of true positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
       FPR: Nx4-elements array containing the rate  of false positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
       PPV: Nx4-elements array containing the Positive Predictive Values that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
       NFP: Number of expected false positives (FPR times the ~7 million stars in the Besancon simulation of the Solar neighborhood) 
       
Each component of the 4-elements dimension of TPR, FPR and PPV corresponds to a different mode of input data, see the description of "LN_PRIOR" above for more detail.

When this fits file is used, the Bayesian probabilities of each star will be associated with a TPR, FPR, NFP and PPV values in the METRICS sub-structure of the output structure.
           
This file must be located at /data/banyan_sigma_metrics.fits in the directory where BANYAN_SIGMA.pro is compiled. The file provided with this release corresponds to the set of models described in Gagné et al. (2018).
           
## CALLING SEQUENCE:
 
```python
OUTPUT_STRUCTURE = BANYAN_SIGMA(stars_data=None, column_names=None, hypotheses=None, ln_priors=None, ntargets_max=None, ra=None, dec=None, pmra=None, pmdec=None, epmra=None, epmdec=None, dist=None, edist=None, rv=None, erv=None, psira=None, psidec=None, epsira=None, epsidec=None, plx=None, eplx=None, constraint_dist_per_hyp=None, constraint_edist_per_hyp=None, unit_priors=True/False, lnp_only=True/False, no_xyz=True/False, use_rv=True/False, use_dist=True/False, use_plx=True/False, use_psi=True/False)
```

## OPTIONAL INPUTS:

       stars_data: An structured array that contains at least the following informations: ra, dec, pmra, pmdec, epmra, and epmdec. It can also optionally contain the informations on rv, erv, dist, edist, plx, eplx, psira, psidec, epsira, epsidec. See the corresponding keyword descriptions for more information. If this input is not used, the keywords ra, dec, pmra, pmdec, epmra, and epmdec must all be specified.
       column_names: A Python dictionary that contains the names of the "stars_data" columns columns which differ from the default values listed above. For example, column_names = {'RA':'ICRS_RA'} can be used to specify that the RA values are listed in the column of stars_data named ICRS_RA.
       ra: Right ascension (decimal degrees). A N-elements array can be specified to calculate the Bayesian probability of several stars at once, but then all mandatory inputs must also be N-elements arrays.
       dec: Declination (decimal degrees).
       pmra: Proper motion in the right ascension direction (mas/yr, must include the cos(dec) factor).
       pmdec: Proper motion in the declination direction (mas/yr).
       epmra: Measurement error on the proper motion in the right ascension direction (mas/yr, must not include the cos(dec) factor).
       epmdec:  Measurement error on the proper motion in the declination direction (mas/yr).
       rv: Radial velocity measurement to be included in the Bayesian probability (km/s). If this keyword is set, erv must also be set. A N-elements array must be used if N stars are analyzed at once.
       erv: Measurement error on the radial velocity to be included in the Bayesian probability (km/s). A N-elements array must be used if N stars are analyzed at once.
       dist: Distance measurement to be included in the Bayesian probability (pc). By default, the banyan_sigma() Bayesian priors are meant for this keyword to be used with trigonometric distances only. Otherwise, the rate of true positives may be far from the nominal values described in Gagné et al. (2018). If this keyword is set, edist must also be set. A N-elements array must be used if N stars are analyzed at once.
       edist: Measurement error on the distance to be included in the Bayesian probability (pc). A N-elements array must be used if N stars are analyzed at once.
       plx: Parallax measurement to be included in the Bayesian probability (mas). The distance will be approximated with dist = 1000/plx. If this keyword is set, eplx must also be set. A N-elements array must be used if N stars are analyzed at once.
       eplx: Measurement error on the parallax to be included in the Bayesian probability (mas). The distance error will be approximated with edist = 1000/plx**2*eplx. A N-elements array must be used if N stars are analyzed at once.
       psira: Parallax motion factor PSIRA described in Gagné et al. (2018), in units of 1/yr. If this keyword is set, the corresponding psidec, epsira and epsidec keywords must also be set. This measurement is only useful when proper motions are estimated from two single-epoch astrometric measurements. It captures the dependence of parallax motion as a function of distance, and allows banyan_sigma() to shift the UVW center of the moving group models, which is equivalent to correctly treating the input "proper motion" pmra, pmdec, epmra, epmdec as a true apparent motion. This keyword should *not* be used if proper motions were derived from more than two epochs, or if they were obtained from a full parallax solution. A N-elements array must be used if N stars are analyzed at once.
       psidec: Parallax motion factor psidec described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
       epsira: Measurement error on the parallax motion factor psira described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
       epsidec: Measurement error on the parallax motion factor psidec described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
       ntargets_max: (default 10^6). Maximum number of objects to run at once in BANYAN_SIGMA to avoid saturating the RAM. If more targets are supplied, banyan_sigma() runs over a loop of several batches of ntargets_max objects.
       hypotheses: The list of Bayesian hypotheses to be considered. They must all be present in the parameters fits file (See REQUIREMENTS #1 above).
       ln_priors: An dictionary that contains the natural logarithm of Bayesian priors that should be *multiplied with the default priors* (use unit_priors=True if you want only ln_priors to be considered). The structure must contain the name of each hypothesis as keys, and the associated scalar value of the natural logarithm of the Bayesian prior for each key. 
       constraint_dist_per_hyp: A structured array that contains a distance constraint (in pc). Each of the Bayesian hypotheses must be included as keys and the distance must be specified as its associated scalar value. constraint_edist_per_hyp must also be specified if constraint_dist_per_hyp is specified. This keyword is useful for including spectro-photometric distance constraints that depend on the age of the young association or field.
       constraint_edist_per_hyp: A structured array that contains a measurement error on the distance constraint (in pc). Each of the Bayesian hypotheses must be included as keys and the distance error must be specified as its associated scalar value.

## OPTIONAL INPUT KEYWORDS:

       unit_priors: If this keyword is set, all default priors are set to 1 (but they are still overrided by manual priors input with the keyword ln_priors).
       lnp_only: If this keyword is set, only Bayesian probabilities will be calculated and returned.
       no_xyz: If this keyword is set, the width of the spatial components of the multivariate Gaussian will be widened by a large factor, so that the XYZ components are effectively ignored. This keyword must be used with extreme caution as it will generate a significant number of false-positives and confusion between the young associations.
       use_rv: Use any radial velocity values found in the stars_data input structure.
       use_dist: Use any distance values found in the stars_data input structure.
       use_plx: Use any parallax values found in the stars_data input structure.
       use_psi: Use any psira, psidec values found in the stars_data input structure.

## OUTPUT:

This routine outputs a structured array, with the following keys:

       NAME: The name of the object (as taken from the input structure).
       ALL: A structure that contains the Bayesian probability (0 to 1) for each of the associations (as individual keys).
       METRICS: A structure that contains the performance metrics associated with the global Bayesian probability of this target. This sub-structure contains the following keys:
              TPR: Rate of true positives expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
              FPR: Rate of false positives (from the field) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
              PPV: Positive Predictive Value (sample contamination) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
       BESTYA_STR: A sub-structure similar to those described above for the most probable young association (ignoring the field possibility).
       YA_PROB: The Bayesian probability (0 to 1) that this object belongs to any young association (i.e., excluding the field).
       LIST_PROB_YAS: A list of young associations with at least 5% Bayesian probability. Their relative probabilities (%) are specified between parentheses.
       BEST_HYP: Most probable Bayesian hypothesis (including the field)
       BEST_YA: Most probable single young association.
       [ASSOCIATION_1]: Sub-structure containing the relevant details for assiciation [ASSOCIATION_1].
       [ASSOCIATION_2]: (...)
       (...)
       [ASSOCIATION_N] - (...)

These per-association sub-structures contain the following keys:

       HYPOTHESIS: Name of the association.
       LN_P: Natural logarithm of the Bayesian probability (LN of 0 to 1).
       D_OPT: Optimal distance (pc) that maximizes the Bayesian likelihood for this hypothesis.
       RV_OPT: Optimal radial velocity (km/s) that maximizes the Bayesian likelihood for this hypothesis.
       ED_OPT: Error on the optimal distance (pc), which approximates the 68% width of how the likelihood varies with distance.
       ERV_OPT: Error on the optimal radial velocity (km/s), which approximates the 68% width of how the likelihood varies with radial velocity.
       XYZUVW: 6-dimensional array containing the XYZ and UVW position of the star at the measured radial velocity and/or distance, or the optimal radial velocity and/or distance when the first are not available (units of pc and km/s).
       EXYZUVW: Errors on XYZUVW (units of pc and km/s).
       XYZ_SEP: Separation between the optimal or measured XYZ position of the star and the center of the multivariate Gaussian model of this Bayesian hypothesis (pc).
       UVW_SEP: Separation between the optimal or measured UVW position of the star and the center of the multivariate Gaussian model of this Bayesian hypothesis (km/s).
       XYZ_SEP: N-sigma separation between the optimal or measured XYZ position of the star and the multivariate Gaussian model of this Bayesian hypothesis (no units).
       UVW_SEP: N-sigma separation between the optimal or measured UVW position of the star and the multivariate Gaussian model of this Bayesian hypothesis (no units).
       MAHALANOBIS: Mahalanobis distance between the optimal or measured XYZUVW position of the star and the multivariate Gaussian model. A Mahalanobis distance is a generalization of a 6D N-sigma distance that accounts for covariances. 

### KNOWN BUGS
* Crashes when no metrics.fits file is provided, but that file should be optional.
* NO_XYZ has not been implemented yet.
