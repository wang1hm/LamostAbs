#IDL to python
# from qsofitmore import QSOFitNew
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits

from scipy.interpolate import interp1d as sp_interp1d
from scipy.optimize import curve_fit

class Filter():
    """A class for the filter's property:"""
    def __init__(self, effec_wave, lam, RC):
        """Initialize the Filter object
        
        Parameters
        ----------
        effec_wave : 1-D np.array
            effective wave of filters
        lam : 3-D np.array
            lam_bin in filters
        RC : 3-D np.array
            responding constant of filters
        """
        self.effec_wave = self.effec_wave
        self.wave = self.wave
        self.RC = self.RC




class CaliFlux():
    """A class for calibration of the flux"""
    def __init__(self, wave, flux, err, mag, emag, z=None, ra=None, dec=None, plateid=None, 
                 mjd=None, fiberid=None, and_mask=None, or_mask=None):
        """Initialize the CaliFlux object.

        Parameters
        ----------
        wave : 1-D np.array
            Wavelength axis of the spectrum
        flux : 1-D np.array
            (Relative) flux density of the spectrum
        err : 1-D np.array
            (Relative) flux error (one sigma) of the spectrum
        mag : 1-D np.array
            magnitude in g,r,i band
        emag : 1-D np.array
            error of magnitude in g,r,i band
        z : float, optional
            Redshift, by default None
        ra : float, optional
            Right ascension of the object, by default None
        dec : float, optional
            Declination of the object, by default None
        plateid : int, optional
            LAMOST plateid of the object, by default None
        mjd : int, optional
            MJD of the observation, by default None
        fiberid : int, optional
            LAMOST fiberid of the object, by default None
        and_mask : 1-D np.array, optional
            And mask of the flux, by default None
        or_mask : 1-D np.array, optional
            Or mask of the flux, by default None
        """
        self.wave = np.asarray(wave, dtype=np.float64)
        self.flux = np.asarray(flux, dtype=np.float64)
        self.err = np.asarray(err, dtype=np.float64)
        self.sn_obs = self.flux/self.err
        self.mag = self.mag
        self.emag = self.emag
        self.z = z
        self.and_mask = and_mask
        self.or_mask = or_mask
        self.ra = ra
        self.dec = dec
        # self.name = name
        self.plateid = plateid
        self.mjd = mjd
        self.fiberid = fiberid

    def mag2flux(mag,emag):
        """transform the magnitude into flux"""
        effec_wave = Filter.effec_wave
        flam = 10**(-0.4 * (mag + 2.406 + 5.*np.log10(effec_wave)))
        eflam = flam * np.log(10.) * 0.4 * emag
        return flam,eflam

    def fit_p_RedBlue(flux_fit,p0,p1):
        """the fit function to optimize the calibration - RedBlue mode"""
        res = np.zeros(3)
        res[0]   = flux_fit[0]*p0
        res[1:3] = np.array(flux_fit[1:3])*p1
        return res
        
    def fit_p_Together(flux_fit,p):
        """the fit function to optimize the calibration - RedBlue mode"""
        res = np.zeros(3)
        res = flux_fit * p
        return res
    
    