#IDL to python
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
#from QSOFit import QSOFit
COLOR = ['#5f160d','#d85c47','#cf4733','#f26c28','#f7986c','#f3b63e',
       '#ffbb00','#fde24a','#fff988','#c5e127','#affd4a','#74d012',
       '#d4ff78','#0baf46','#99f6d1','#30b494','#9ffcff','#22cdeb',
       '#2bffe9','#7298c5','#176ab1','#2348e4','#888ad8','#3109cc',
       '#b873f1','#8224f3','#e0aaf2','#f224f3','#eb89a8','#ef3762',
       '#696969','#ff0030']

#Base on IDL code
""""
function cmd,p
   common LAMOSTspec,SumFlux_g,SumFlux_r,SumFlux_i,SumFlux_z;,sdss_effec_wave
   ;sdss_effec_wave = [4686.,6165.,7481.,8931.]
   flux_fit = dblarr(3)
   flux_fit[0]   = SumFlux_g*p(0)
   flux_fit[1:2] =[SumFlux_r,SumFlux_i]*p(1)
   print,flux_fit
   return,flux_fit
end
"""

#sdss_effec_wave
float LAMOSTspec,SumFlux_g,SumFlux_r,SumFlux_i,SumFlux_z

def cmd(p):
    #sdss_effec_wave = [4686.,6165.,7481.,8931.]
    flux_fit = np.zeros(3)
    flux_fit[0]   = SumFlux_g*p[0]
    flux_fit[1:2] =[SumFlux_r,SumFlux_i]*p[1]
    print(flux_fit)
    return flux_fit

"""
;---------------read Filter-------------------
readcol,'test.dat',basename,obsid,RA,DEC,psfMag_g,psfMag_r,psfMag_i,psfMag_z,psfMagErr_g,psfMagErr_r,psfMagErr_i,psfMagErr_z,together,Photometry,redshift,f='(a,d,d,d,d,d,d,d,d,d,d,d,x,a,a,d)'

num             = n_elements(basename)

"""

data = pd.read_csv('./dr45_match_file.csv')
num = len(data)

fig, ax = plt.subplots(2,1)

Emission_name = ['Ha','Hb','MgII','[OIII]',' ','[NII]',' ','[SII]',' ','CIII','CIV','Lya']
Emission_wave = [6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216]*(1.+redshift[i])

SDSS_path = './Filter/SDSS/'
mag_sdss_g = pd.read_csv(SDSS_path+'g.dat')
mag_sdss_r = pd.read_csv(SDSS_path+'r.dat')
mag_sdss_i = pd.read_csv(SDSS_path+'i.dat')
mag_sdss_z = pd.read_csv(SDSS_path+'z.dat')
mag_sdss_g.iloc=lam_g*10.
lam_r=lam_r*10.
lam_i=lam_i*10.
lam_z=lam_z*10.

Pstr_path = './Filter/Panstarr/'
mag_pstr_g = pd.read_csv(Pstr_path+'g.dat')
mag_pstr_r = pd.read_csv(Pstr_path+'r.dat')
mag_pstr_i = pd.read_csv(Pstr_path+'i.dat')
mag_pstr_z = pd.read_csv(Pstr_path+'z.dat')

def RedBlue(p):
    
