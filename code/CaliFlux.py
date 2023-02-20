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

#sdss_effec_wave
float LAMOSTspec,SumFlux_g,SumFlux_r,SumFlux_i,SumFlux_z

def cmd(p):
    #sdss_effec_wave = [4686.,6165.,7481.,8931.]
    flux_fit = np.zeros(3)
    flux_fit[0]   = SumFlux_g*p[0]
    flux_fit[1:2] =[SumFlux_r,SumFlux_i]*p[1]
    print(flux_fit)
    return flux_fit

data = pd.read_csv('./dr45_match_file.csv')
num = len(data)

fig, ax = plt.subplots(2,1)

Emission_name = ['Ha','Hb','MgII','[OIII]',' ','[NII]',' ','[SII]',' ','CIII','CIV','Lya']
#Emission_wave = [6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216]*(1.+redshift[i])
def read_mag(photometry):
    if photometry == 'SDSS':
        SDSS_path = './Huimei/Filter/SDSS/'
        mag_sdss_g = pd.read_csv(SDSS_path+'g.dat',names=['Del','lam_g','RC_g'],sep=' ')
        lam_g = mag_sdss_g['lam_g']
        mag_sdss_r = pd.read_csv(SDSS_path+'r.dat',names=['Del','lam_r','RC_r'],sep=' ')
        lam_r = mag_sdss_r['lam_r']
        mag_sdss_i = pd.read_csv(SDSS_path+'i.dat',names=['Del','lam_i','RC_i'],sep=' ')
        lam_i = mag_sdss_i['lam_i']
        mag_sdss_z = pd.read_csv(SDSS_path+'z.dat',names=['Del','lam_z','RC_z'],sep=' ')
        lam_z = mag_sdss_z['lam_z']
        lam_g=lam_g*10.
        lam_r=lam_r*10.
        lam_i=lam_i*10.
        lam_z=lam_z*10.
    if photometry == 'Panstarr':
        PSTR_path = './Huimei/Filter/Panstarr/'
        mag_pstr_g = pd.read_csv(PSTR_path+'g.dat',names=['Del','lam_g','RC_g'],sep=' ')
        lam_g = mag_ptsr_g['lam_g']
        mag_pstr_r = pd.read_csv(PSTR_path+'r.dat',names=['Del','lam_r','RC_r'],sep=' ')
        lam_r = mag_ptsr_r['lam_r']
        mag_pstr_i = pd.read_csv(PSTR_path+'i.dat',names=['Del','lam_i','RC_i'],sep=' ')
        lam_i = mag_ptsr_i['lam_i']
        mag_pstr_z = pd.read_csv(PSTR_path+'z.dat',names=['Del','lam_z','RC_z'],sep=' ')
        lam_z = mag_ptsr_z['lam_z']
        lam_g=lam_g*10.
        lam_r=lam_r*10.
        lam_i=lam_i*10.
        lam_z=lam_z*10.
    return lam_g, lam_r, lam_i, lam_z

def RedBlue(p,z,photometry,name):
    Emission_wave = np.array([6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216])*(1.+z)
    specdata = fits.open('Volumes/My_Passport/Calibration/lamost_dr45/'+name+'.fits.gz')
    lam_g, lam_r, lam_i, lam_z = read_mag(photometry)
    flux = specdata[0].data[2]#1e-17
    ivar = specdata[0].data[1]
    lam = specdata[0].data[2]
    andmask = specdata[0].data[3]
    ormask = specdata[0].data[4]
    err = 1./np.sqrt(ivar)
    wave_index1 = np.where((3900<lam) &(lam < 4100))
    wave_index2 = np.where((7900<lam) &(lam < 8100))
    flux1        = median(flux[wave_index1])
    flux2        = median(flux[wave_index2])
    flux_max     = max([flux1,flux2])
    flux_min     = min([flux1,flux2])
    lam_min=int(min(lam))+1
    lam_max=int(max(lam))
    lam_bin=np.arange(0,lam_max-lam_min,1)+lam_min
    #set the bin for griz band
    lam_bin_g = lam_bin[np.where((min(lam_g)<lam_bin)&(lam_bin < max(lam_g)))]
    lam_bin_r = lam_bin[np.where((min(lam_r)<lam_bin)&(lam_bin < max(lam_r)))]
    lam_bin_i = lam_bin[np.where((min(lam_i)<lam_bin)&(lam_bin < max(lam_i)))]
    lam_bin_z = lam_bin[np.where((min(lam_z)<lam_bin)&(lam_bin < max(lam_z)))]

