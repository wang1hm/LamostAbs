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
        self.effec_wave = effec_wave
        self.lam = lam
        self.RC = RC


class CaliFlux():
    """A class for calibration of the flux"""
    def __init__(self, lam, flux, err, mag, emag, z=None, ra=None, dec=None, plateid=None, 
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
        self.lam = np.asarray(lam, dtype=np.float64)
        self.flux = np.asarray(flux, dtype=np.float64)
        self.err = np.asarray(err, dtype=np.float64)
        self.sn_obs = self.flux/self.err
        self.mag = mag
        self.emag = emag
        self.z = z
        self.and_mask = and_mask
        self.or_mask = or_mask
        self.ra = ra
        self.dec = dec
        # self.name = name
        self.plateid = plateid
        self.mjd = mjd
        self.fiberid = fiberid

    def mag2flux(self):
        """transform the magnitude into flux"""
        effec_wave = Filter.effec_wave
        flam = 10**(-0.4 * (self.mag + 2.406 + 5.*np.log10(effec_wave)))
        eflam = flam * np.log(10.) * 0.4 * self.emag
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
    
    def RedBlue(self):
        """the calibration code for RedBlue mode"""
        psfMag_g,psfMag_r,psfMag_i,psfMag_z = self.mag[0],self.mag[1],self.mag[2],self.mag[3]
        psfMagErr_g, psfMagErr_r, psfMagErr_i, psfMagErr_z = self.emag[0],self.emag[1],self.emag[2],self.emag[3]
        lam_g, lam_r, lam_i, lam_z, RC_g, RC_r, RC_i, RC_z = Filter.lam[0],Filter.lam[1],Filter.lam[2],Filter.mag[3],Filter.RC[0],Filter.RC[1],Filter.RC[2],Filter.Rc[3]
        wave_index1 = np.where((3900<self.lam) &(self.lam < 4100))
        wave_index2 = np.where((7900<self.lam) &(self.lam < 8100))
        flux1 = np.median(self.flux[wave_index1])
        flux2 = np.median(self.flux[wave_index2])
        flux_max = np.max([flux1,flux2])
        flux_min = np.min([flux1,flux2])
        lam_min=int(min(self.lam))+1
        lam_max=int(max(self.lam))
        lam_bin=np.arange(0,lam_max-lam_min,1)+lam_min

        #set the bin for griz band
        lam_bin_g = lam_bin[np.where((min(lam_g)<lam_bin)&(lam_bin < max(lam_g)))]
        lam_bin_r = lam_bin[np.where((min(lam_r)<lam_bin)&(lam_bin < max(lam_r)))]
        lam_bin_i = lam_bin[np.where((min(lam_i)<lam_bin)&(lam_bin < max(lam_i)))]
        #lam_bin_z = lam_bin[np.where((min(lam_z)<lam_bin)&(lam_bin < max(lam_z)))]

        #set the bin for RC
        f_g = sp_interp1d(lam_g, RC_g)
        RC_g_bin = f_g(lam_bin_g)
        f_r = sp_interp1d(lam_r, RC_r)
        RC_r_bin = f_r(lam_bin_r)
        f_i = sp_interp1d(lam_i, RC_i)
        RC_i_bin = f_i(lam_bin_i)
        #f_z = sp_interp1d(lam_z, RC_z)
        #RC_z_bin = f_z(lam_bin_z)

        #interpolate to get the flux for lambda bin
        f_fl = sp_interp1d(self.lam,self.flux)
        flux_bin_g = f_fl(lam_bin_g)
        flux_bin_r = f_fl(lam_bin_r)
        flux_bin_i = f_fl(lam_bin_i)
        #flux_bin_z = f_fl(lam_bin_z)

        SumFlux_g   = np.sum(RC_g_bin*flux_bin_g*lam_bin_g)/np.sum(RC_g_bin*lam_bin_g)
        SumFlux_r   = np.sum(RC_r_bin*flux_bin_r*lam_bin_r)/np.sum(RC_r_bin*lam_bin_r)
        SumFlux_i   = np.sum(RC_i_bin*flux_bin_i*lam_bin_i)/np.sum(RC_i_bin*lam_bin_i)
        #SumFlux_z   = np.sum(RC_z_bin*flux_bin_z*lam_bin_z)/np.sum(RC_z_bin*lam_bin_z)


