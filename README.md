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

The latest version of this package was only tested with Python 3.12.3 installed with a Mac M1-native conda setup.

We highly recommend using a Python virtual environment (or a conda environment) to use this package, as some Python packages may have issued non-retro compatible updates since the publication of this package. This can be done with the following steps once you are in the banyan_sigma directory (or anywhere else you prefer).

       python -m venv banyan_sigma_env

This will initiate a virtual environment where your Python packages specific to banyan_sigma will live, without disrupting your other installations.

Once this is done, you need to activate this virtual environment with the following terminal command:

       source banyan_sigma_env/bin/activate
	   
You can then install banyan sigma with the following command:

       pip install git+https://github.com/jgagneastro/banyan_sigma.git

This should now install the exact versions of all the required subpackages. However, if you have locally downloaded the BANYAN Sigma tool and included it in your Python Path environment variable, you may need to do this manually, using:

       pip install -r [your/own/path]/banyan_sigma/requirements.txt
       
If you prefer to use conda instead of pyenv, you can follow these steps:

       conda create --name bsigma_env python==3.12.3
       conda activate bsigma_env
       pip install git+https://github.com/jgagneastro/banyan_sigma.git

Note that you can choose your own environment name instead of bsigma_env.

## A SIMPLE EXAMPLE:

Once the packages are installed, you can start using banyan_sigma. First open Python:

       python

Here is an example code to calculate a membership probability using a specific set of observables (for the star AU Mic) which you can paste directly in Python:
       
    from banyan_sigma import *
    
    #Define observables for an example star
    ra=311.2911826481039
    dec=-31.3425000799281
    
    #Proper motions are provided in mas/yr, pmra is implicitly pmra*cos(dec)
    pmra=281.319
    epmra=0.022
    pmdec=-360.148
    epmdec=0.019
    
    #Parallaxes are provided in mas
    plx=102.943
    eplx=0.023
    
    #Radial velocities are provided in km/s
    rv=-5.2
    erv=0.7
    
    #Determine membership probability
    output = membership_probability(ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,epmra=epmra,epmdec=epmdec,plx=plx,eplx=eplx,rv=rv,erv=erv, use_plx=True, use_rv=True)
    
    #Inverstigate the outputs
    output.iloc[0]

Because the output is a Pandas dataframe with 1 row (as we calculated the membership probabilities for one star), one might output the overall results with:

    >>> output.iloc[0]
    118TAU         LN_P      -74.990764
                   D_OPT       9.714114
                   RV_OPT          -5.2
                   ED_OPT       0.00217
                   ERV_OPT          0.7
                                ...    
    ALL            FIELD       0.001812
    YA_PROB        Global      0.998188
    LIST_PROB_YAS  Global          BPMG
    BEST_HYP       Global          BPMG
    BEST_YA        Global          BPMG
    Name: 0, Length: 763, dtype: object

Detailed membership probabilities (values go from zero to one) for all associations may be output with the following:

```python
output['ALL'].iloc[0]

118TAU     1.016362e-29
ABDMG      9.260548e-40
BPMG       9.981879e-01
CAR        2.304203e-27
CARN       1.120327e-14
CBER       0.000000e+00
COL        1.317175e-31
CRA        0.000000e+00
EPSC      8.142372e-234
ETAC       0.000000e+00
HYA        1.365976e-95
IC2391    1.093534e-188
IC2602     0.000000e+00
LCC        4.385118e-17
OCT       3.659492e-183
PL8       1.310888e-273
PLE        2.264604e-72
ROPH       0.000000e+00
TAU        1.364748e-34
THA        1.722218e-78
THOR       0.000000e+00
TWA        8.462222e-24
UCL        8.794057e-12
UCRA      3.927025e-237
UMA        0.000000e+00
USCO       1.383815e-20
XFOR      8.764554e-118
VCA       1.573234e-102
ARG        8.089683e-20
PRAE       0.000000e+00
MUTAU      1.241229e-21
PERI       3.295221e-75
FIELD      1.812123e-03
Name: 0, dtype: float64
```

Detailed properties for a given association such as Beta Pictoris (BPMG) can be obtained with:

```python
output['BPMG'].iloc[0]

LN_P           -8.233839
D_OPT           9.714114
RV_OPT         -5.200000
ED_OPT          0.002170
ERV_OPT         0.700000
X               7.589092
Y               1.703744
Z              -5.819531
U             -10.563371
V             -16.192239
W              -9.835885
EX              0.001696
EY              0.000381
EZ              0.001300
EU              0.546874
EV              0.122828
EW              0.419364
XYZ_SEP        13.427037
UVW_SEP         0.886663
XYZ_SIG         1.721692
UVW_SIG         0.874458
MAHALANOBIS     1.748404
Name: 0, dtype: float64

```

(2) A fits file containing the parameters of the multivariate Gaussian models of each Bayesian hypothesis must be included at /data/banyan_sigma_parameters.fits in the directory where banyan_sigma_ is compiled. The file provided with this release corresponds to the set of 27 young associations described in Gagné et al. (2018). The fits file can be written with the IDL MWRFITS.PRO function from an IDL array of structures of N elements, where N is the total number of multivariate Gaussians used in the models of all Bayesian hypotheses. Each element of this structure contains the following information:

       NAME: The name of the model (scalar string).
       CENTER_VEC: Central XYZUVW position of the model (6D vector, in units of pc and km/s).
       COVARIANCE_MATRIX: Covariance matrix in XYZUVW associated with the model (6x6 matrix, in mixed units of pc and km/s).
       PRECISION_MATRIX: (Optional) Matrix inverse of COVARIANCE_MATRIX, to avoid re-calculating it many times (6x6 matrix).
       LN_NOBJ: (Optional) Natural logarithm of the number of objects used to build the synthetic model (scalar). This is not used in banyan_sigma.
       COVARIANCE_DETERM: (Optional) Determinant of the covariance matrix, to avoid re-calculating it many times (scalar).
       PRECISION_DETERM: (Optional) Determinant of the precision matrix, to avoid re-calculating it many times (scalar).
       LN_ALPHA_K: (Optional) Natural logarithm of the alpha_k inflation factors that ensured a fixed rate of true positives at a given Bayesian probability treshold. See Gagné et al. (2018) for more detail (scalar or 4-elements vector). This is not used in BANYAN_SIGMA.
       LN_PRIOR: Natural logarithm of the Bayesian prior (scalar of 4-elements vector). When this is a 4-elements vector, the cases with only proper motion, proper motion + radial velocity, proper motion + distance or proper motion + radial velocity + distance will be used with the corresponding element of the LN_PRIOR vector.
       LN_PRIOR_OBSERVABLES: Scalar string or 4-elements vector describing the observing modes used for each element of ln_prior. This is not used in banyan_sigma.
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

Currently, we have not provided this file because it is rarely used, and the machine-learning performance metrics have not been computed for some of the more recent young associations.
           
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
