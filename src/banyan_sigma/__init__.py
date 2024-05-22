__author__ = 'Jonathan Gagne'
__email__ = 'jonathan.gagne@astro.umontreal.ca'
__uri__ = "https://github.com/jgagneastro/banyan_sigma"
__license__ = "MIT"
__description__ = "A Bayesian classifier for young moving groups membership"

#Moca should be available automatically when importing mocapy
__all__ = ["membership_probability"]

#Import the local moca.py file
from .core import membership_probability