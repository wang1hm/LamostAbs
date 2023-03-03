#!/usr/bin/env python
from astropy.io import fits
import numpy as np
import pandas as pd
from glob import glob
import re
import matplotlib.pyplot as plt
import os
# from PyAstronomy import pyasl
from scipy.signal import savgol_filter
from scipy.stats import sigmaclip
from pathlib import Path
from astropy.nddata import StdDevUncertainty,VarianceUncertainty,InverseVariance
from astropy.table import Table
from astropy import units as u
# from specutils import Spectrum1D,SpectrumCollection,SpectrumList
# from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler, median_smooth
import warnings
from astropy.units import Quantity

class LowResSpec():
    def __init__(self, wave, flux, err, z=None, ra=None, dec=None, plateid=None, 
                 mjd=None, fiberid=None, and_mask=None, or_mask=None):
        self.wave = np.asarray(wave, dtype=np.float64)
        self.flux = np.asarray(flux, dtype=np.float64)
        self.err = np.asarray(err, dtype=np.float64)
        self.sn_obs = self.flux/self.err
        self.z = z
        self.and_mask = and_mask
        self.or_mask = or_mask
        self.ra = ra
        self.dec = dec
        # self.name = name
        self.plateid = plateid
        self.mjd = mjd
        self.fiberid = fiberid


    @classmethod
    def load(cls, fname):
        hdu = fits.open(fname)
        header = hdu[0].header
        ra=header['RA']
        dec=header['DEC']
        plateid = header['OBSID']
        mjd = header['MJD'] 
        fiberid = header['FIBERID']
        objid = header['OBJNAME']
        objname = header['DESIG']
        data = hdu[0].data
        flux = data[0]
        ivar = data[1]
        wave = data[2]
        and_mask = data[3]
        or_mask = data[4] 
        z_pipe = header['z']
        ivar = pd.Series(data[1])
        ivar.replace(0, np.nan, inplace=True)
        ivar_safe = ivar.interpolate()
        err = 1./np.sqrt(ivar_safe.values) 
        hdu.close()  
        return cls(wave=wave, flux=flux, err=err, z=z_pipe, ra=ra, dec=dec, plateid=plateid, 
                   mjd=mjd, fiberid=fiberid, and_mask=and_mask, or_mask=or_mask)