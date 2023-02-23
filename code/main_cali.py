#Calibrate the Lamost spectra
from qsofitmore import QSOFitNew
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import glob
import os
import shutil
import wget
import ftplib
from scipy.interpolate import interp1d as sp_interp1d
from scipy.optimize import curve_fit
import CaliFlux

#pd.to_numeric(data.iloc[:,5])
mag = data.iloc[0][5:9]
mag_err = data.iloc[0][9:13]
#p = data['redshift'][0]
z = float(data['redshift'][0])
photometry = 'SDSS'
name = data['basename'][0]
Res = RedBlue(mag,mag_err,z,photometry,str(name))
plt.plot(Res[2],Res[0],zorder = 0)
plt.scatter([4686.,6165.,7481.],arr_f1,c='r',zorder=1)