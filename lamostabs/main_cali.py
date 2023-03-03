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
import CaliFlux

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d as sp_interp1d
#from QSOFit import QSOFit
COLOR = ['#5f160d','#d85c47','#cf4733','#f26c28','#f7986c','#f3b63e',
       '#ffbb00','#fde24a','#fff988','#c5e127','#affd4a','#74d012',
       '#d4ff78','#0baf46','#99f6d1','#30b494','#9ffcff','#22cdeb',
       '#2bffe9','#7298c5','#176ab1','#2348e4','#888ad8','#3109cc',
       '#b873f1','#8224f3','#e0aaf2','#f224f3','#eb89a8','#ef3762',
       '#696969','#ff0030']

data = pd.read_csv('./02_Calibration/dr45_match_file.csv',names=['basename','obsid','RA','DEC','psfMag_g',\
        'psfMag_r','psfMag_i','psfMag_z','psfMagErr_g','psfMagErr_r','psfMagErr_i','psfMagErr_z',\
        'together','Photometry','redshift','f'])
data=data.drop(0)
data=data.reset_index(drop=False)
#num = len(data)-lam_min
for i in [5,6,7,8,9,10,11,12,15]:
    data.iloc[:,i]=pd.to_numeric(data.iloc[:,i])

#pd.to_numeric(data.iloc[:,5])
mag = data.iloc[0][5:9]
mag_err = data.iloc[0][9:13]
#p = data['redshift'][0]
z = data['redshift'][0]
photometry = data['Photometry'][0]
name = data['basename'][0]
Res,arr_f1 = CaliFlux.RedBlue(mag,mag_err,z,photometry,str(name))
plt.plot(Res[2],Res[0],zorder = 0)
plt.scatter([4686.,6165.,7481.],arr_f1,c='r',zorder=1)