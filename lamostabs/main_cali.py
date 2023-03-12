# from qsofitmore import QSOFitNew
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import glob
import os
import shutil
import ftplib
import CaliFlux

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d as sp_interp1d
#from QSOFit import QSOFit

data = pd.read_csv('./02_Calibration/dr45_match_file.csv',names=['basename','obsid','RA','DEC','psfMag_g',\
        'psfMag_r','psfMag_i','psfMag_z','psfMagErr_g','psfMagErr_r','psfMagErr_i','psfMagErr_z',\
        'together','Photometry','redshift','f'])
data=data.drop(0)
data=data.reset_index(drop=False)
#num = len(data)-lam_min
for i in [5,6,7,8,9,10,11,12,15]:
    data.iloc[:,i]=pd.to_numeric(data.iloc[:,i])

#pd.to_numeric(data.iloc[:,5])
spec_index=0
mag = data.iloc[spec_index][5:9]
mag_err = data.iloc[spec_index][9:13]
#p = data['redshift'][spec_index]
z = data['redshift'][spec_index]
photometry = data['Photometry'][spec_index]
name = data['basename'][spec_index]
Res,arr_f1 = CaliFlux.RedBlue(mag,mag_err,z,photometry,str(name))


plt.plot(Res[2],Res[0],zorder = 0)
plt.scatter([4686.,6165.,7481.],arr_f1,c='r',zorder=1)

