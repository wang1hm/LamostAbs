#IDL to python
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

from scipy.interpolate import interp1d as sp_interp1d
from scipy.optimize import curve_fit
#from QSOFit import QSOFit
COLOR = ['#5f160d','#d85c47','#cf4733','#f26c28','#f7986c','#f3b63e',
       '#ffbb00','#fde24a','#fff988','#c5e127','#affd4a','#74d012',
       '#d4ff78','#0baf46','#99f6d1','#30b494','#9ffcff','#22cdeb',
       '#2bffe9','#7298c5','#176ab1','#2348e4','#888ad8','#3109cc',
       '#b873f1','#8224f3','#e0aaf2','#f224f3','#eb89a8','#ef3762',
       '#696969','#ff0030']

#sdss_effec_wave
#float LAMOSTspec,SumFlux_g,SumFlux_r,SumFlux_i,SumFlux_z


def cmd(flux_fit,p0,p1):
    #sdss_effec_wave = [4686.,6165.,7481.,8931.]
    res = np.zeros(3)
    res[0]   = flux_fit[0]*p0
    res[1:3] = np.array(flux_fit[1:3])*p1
    return res

#Emission_name = ['Ha','Hb','MgII','[OIII]',' ','[NII]',' ','[SII]',' ','CIII','CIV','Lya']
#Emission_wave = [6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216]*(1.+redshift[i])
def read_mag(photometry):
    if photometry == 'SDSS':
        SDSS_path = './Huimei/Filter/SDSS/'
        mag_sdss_g = pd.read_csv(SDSS_path+'g.dat',names=['Del','lam_g','RC_g'],sep=' ')
        lam_g, RC_g = mag_sdss_g['lam_g'],mag_sdss_g['RC_g']
        mag_sdss_r = pd.read_csv(SDSS_path+'r.dat',names=['Del','lam_r','RC_r'],sep=' ')
        lam_r, RC_r = mag_sdss_r['lam_r'],mag_sdss_r['RC_r']
        mag_sdss_i = pd.read_csv(SDSS_path+'i.dat',names=['Del','lam_i','RC_i'],sep=' ')
        lam_i, RC_i = mag_sdss_i['lam_i'],mag_sdss_i['RC_i']
        mag_sdss_z = pd.read_csv(SDSS_path+'z.dat',names=['Del','lam_z','RC_z'],sep=' ')
        lam_z, RC_z = mag_sdss_z['lam_z'],mag_sdss_z['RC_z']
        lam_g=lam_g*10.
        lam_r=lam_r*10.
        lam_i=lam_i*10.
        lam_z=lam_z*10.
    if photometry == 'Panstarr':
        PSTR_path = './Huimei/Filter/Panstarr/'
        mag_pstr_g = pd.read_csv(PSTR_path+'g.dat',names=['Del','lam_g','RC_g'],sep=' ')
        lam_g, RC_g = mag_ptsr_g['lam_g'],mag_ptsr_g['RC_g']
        mag_pstr_r = pd.read_csv(PSTR_path+'r.dat',names=['Del','lam_r','RC_r'],sep=' ')
        lam_r, RC_r = mag_ptsr_r['lam_r'],mag_ptsr_r['RC_r']
        mag_pstr_i = pd.read_csv(PSTR_path+'i.dat',names=['Del','lam_i','RC_i'],sep=' ')
        lam_i, RC_i = mag_ptsr_i['lam_i'],mag_ptsr_i['RC_i']
        mag_pstr_z = pd.read_csv(PSTR_path+'z.dat',names=['Del','lam_z','RC_z'],sep=' ')
        lam_z, RC_z = mag_ptsr_z['lam_z'],mag_ptsr_z['RC_z']
        lam_g=lam_g*10.
        lam_r=lam_r*10.
        lam_i=lam_i*10.
        lam_z=lam_z*10.
    return lam_g, lam_r, lam_i, lam_z, RC_g, RC_r, RC_i, RC_z

def mag2flux(mag,emag,wave):
    flam = 10**(-0.4 * (mag + 2.406 + 5.*np.log10(wave)))
    eflam = flam * np.log(10.) * 0.4 * emag
    return flam,eflam


def RedBlue(mag,mag_err,z,photometry,name):
    #read the data
    specdata = fits.open('/Volumes/My_Passport/Calibration/lamost_dr45/'+name+'.fits.gz')
    Emission_wave = np.array([6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216])*(1.+z)
    
    psfMag_g,psfMag_r,psfMag_i,psfMag_z = mag[0],mag[1],mag[2],mag[3]
    psfMagErr_g, psfMagErr_r, psfMagErr_i, psfMagErr_z = mag_err[0],mag_err[1],mag_err[2],mag_err[3]
    lam_g, lam_r, lam_i, lam_z, RC_g, RC_r, RC_i, RC_z = read_mag(photometry)
    
    flux = specdata[0].data[0]#1e-17
    ivar = specdata[0].data[1]
    lam = specdata[0].data[2]
    andmask = specdata[0].data[3]
    ormask = specdata[0].data[4]
    err = 1./np.sqrt(ivar)
    wave_index1 = np.where((3900<lam) &(lam < 4100))
    wave_index2 = np.where((7900<lam) &(lam < 8100))
    flux1 = np.median(flux[wave_index1])
    flux2 = np.median(flux[wave_index2])
    flux_max     = max([flux1,flux2])
    flux_min     = min([flux1,flux2])
    lam_min=int(min(lam))+1
    lam_max=int(max(lam))
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

    f_fl = sp_interp1d(lam,flux)
    flux_bin_g = f_fl(lam_bin_g)
    flux_bin_r = f_fl(lam_bin_r)
    flux_bin_i = f_fl(lam_bin_i)
    #flux_bin_z = f_fl(lam_bin_z)

    SumFlux_g   = np.sum(RC_g_bin*flux_bin_g*lam_bin_g)/np.sum(RC_g_bin*lam_bin_g)
    SumFlux_r   = np.sum(RC_r_bin*flux_bin_r*lam_bin_r)/np.sum(RC_r_bin*lam_bin_r)
    SumFlux_i   = np.sum(RC_i_bin*flux_bin_i*lam_bin_i)/np.sum(RC_i_bin*lam_bin_i)
    #SumFlux_z   = np.sum(RC_z_bin*flux_bin_z)/np.sum(RC_z_bin)

    #calculate the flux and error from mag to draw the cali_point
    #arr_f1     = np.zeros(4)
    #arr_f1_err = np.zeros(4)
    arr_f1     = np.zeros(3)
    arr_f1_err = np.zeros(3)
    if photometry == 'SDSS':
        effec_wave = [4686.,6165.,7481.,8931.]
        arr_f1[0],arr_f1_err[0] = np.array(mag2flux(psfMag_g+0.012,psfMagErr_g,wave=effec_wave[0]))/1.0e-17
        arr_f1[1],arr_f1_err[1] = np.array(mag2flux(psfMag_r+0.010,psfMagErr_r,wave=effec_wave[1]))/1.0e-17
        arr_f1[2],arr_f1_err[2] = np.array(mag2flux(psfMag_i+0.028,psfMagErr_i,wave=effec_wave[2]))/1.0e-17
    
    if photometry == 'Panstarr':
        effec_wave = [4810.,6170.,7520.,8660.]
        arr_f1[0],arr_f1_err[0] = np.array(mag2flux(psfMag_g,psfMagErr_g,wave=effec_wave[0]))/1.0e-17
        arr_f1[1],arr_f1_err[1] = np.array(mag2flux(psfMag_r,psfMagErr_r,wave=effec_wave[1]))/1.0e-17
        arr_f1[2],arr_f1_err[2] = np.array(mag2flux(psfMag_i,psfMagErr_i,wave=effec_wave[2]))/1.0e-17
    #arr_f1[3],arr_f1_err[3] = np.array(mag2flux(psfMag_z+0.040,psfMagErr_z,wave=effec_wave[3]))/1.0e-17
   	
    #curve_fit
    params, err = curve_fit(cmd,[SumFlux_g,SumFlux_r,SumFlux_i], arr_f1, sigma = arr_f1_err)

    index_blue = np.where(lam < 5900)
    index_red = np.where(lam > 5900)
    lam_red = lam[index_red]
    lam_blue = lam[index_blue]

    flux1 = np.median(flux[wave_index1]*params[0])
    flux2 = np.median(flux[wave_index2]*params[1])
    flux_max = np.max([flux1,flux2])
    flux_min = np.min([flux1,flux2])

    Res_lam = np.hstack((lam_blue,lam_red))
    Res_flux = np.hstack((flux[index_blue]*params[0],flux[index_red]*params[1]))
    #make the S/N in red part same with the blue part
    #Res_ivar = np.hstack((ivar[index_blue]/(params[0]**2.),ivar[index_red]/(params[1]**2.)))
    Res_ivar = np.hstack((ivar[index_blue]/(params[0]**2.),ivar[index_red]/(params[1]**2.)))
    #Res_ivar = np.hstack((ivar[index_blue],ivar[index_red]))
    #Res_ivar = [ivar[index_blue]/(params[0]^2.),ivar[index_red]/(params[1]^2.)]
    index_ivar0 = np.where(Res_ivar == 0)
    #Res_ivar[index_ivar0] = 0.1
    index_ivar1 = np.where((Res_lam >5700)&(Res_lam < 5900))
    Res_ivar[index_ivar1] = 0.01
    Res_andmask = np.hstack((andmask[index_blue],andmask[index_red]))
    Res_ormask  = np.hstack((ormask[index_blue],ormask[index_red]))
    Res_num = len(Res_lam)
    Res_data    = np.zeros([5,Res_num])
    Res_data[0] = Res_flux
    Res_data[1] = Res_ivar
    Res_data[2] = Res_lam
    Res_data[3] = Res_andmask
    Res_data[4] = Res_ormask
    
    return Res_data, arr_f1, lam_blue, lam_red

def cmd_together(flux_fit,p):
    #sdss_effec_wave = [4686.,6165.,7481.,8931.]
    res = np.zeros(3)
    #res[0]   = flux_fit[0]*p0
    #res[1:3] = np.array(flux_fit[1:3])*p1
    res[0:3] = flux_fit * p
    return res

def Together(mag,mag_err,z,photometry,name):
    #read the data
    specdata = fits.open('/Volumes/My_Passport/Calibration/lamost_dr45/'+name+'.fits.gz')
    Emission_wave = np.array([6563.,4861.,2798,4959,5007,6548,6584,6716,6731,1909,1549,1216])*(1.+z)
    
    psfMag_g,psfMag_r,psfMag_i,psfMag_z = mag[0],mag[1],mag[2],mag[3]
    psfMagErr_g, psfMagErr_r, psfMagErr_i, psfMagErr_z = mag_err[0],mag_err[1],mag_err[2],mag_err[3]
    lam_g, lam_r, lam_i, lam_z, RC_g, RC_r, RC_i, RC_z = read_mag(photometry)
    
    flux = specdata[0].data[0]#1e-17
    ivar = specdata[0].data[1]
    lam = specdata[0].data[2]
    andmask = specdata[0].data[3]
    ormask = specdata[0].data[4]
    err = 1./np.sqrt(ivar)
    wave_index1 = np.where((3900<lam) &(lam < 4100))
    wave_index2 = np.where((7900<lam) &(lam < 8100))
    flux1 = np.median(flux[wave_index1])
    flux2 = np.median(flux[wave_index2])
    flux_max     = max([flux1,flux2])
    flux_min     = min([flux1,flux2])
    lam_min=int(min(lam))+1
    lam_max=int(max(lam))
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

    f_fl = sp_interp1d(lam,flux)
    flux_bin_g = f_fl(lam_bin_g)
    flux_bin_r = f_fl(lam_bin_r)
    flux_bin_i = f_fl(lam_bin_i)
    #flux_bin_z = f_fl(lam_bin_z)

    SumFlux_g   = np.sum(RC_g_bin*flux_bin_g*lam_bin_g)/np.sum(RC_g_bin*lam_bin_g)
    SumFlux_r   = np.sum(RC_r_bin*flux_bin_r*lam_bin_r)/np.sum(RC_r_bin*lam_bin_r)
    SumFlux_i   = np.sum(RC_i_bin*flux_bin_i*lam_bin_i)/np.sum(RC_i_bin*lam_bin_i)
    #SumFlux_z   = np.sum(RC_z_bin*flux_bin_z*lam_bin_z)/np.sum(RC_z_bin*lam_bin_z)

    #calculate the flux and error from mag to draw the cali_point
    #arr_f1     = np.zeros(4)
    #arr_f1_err = np.zeros(4)
    arr_f1     = np.zeros(3)
    arr_f1_err = np.zeros(3)
    if photometry == 'SDSS':
        effec_wave = [4686.,6165.,7481.,8931.]
        arr_f1[0],arr_f1_err[0] = np.array(mag2flux(psfMag_g+0.012,psfMagErr_g,wave=effec_wave[0]))/1.0e-17
        arr_f1[1],arr_f1_err[1] = np.array(mag2flux(psfMag_r+0.010,psfMagErr_r,wave=effec_wave[1]))/1.0e-17
        arr_f1[2],arr_f1_err[2] = np.array(mag2flux(psfMag_i+0.028,psfMagErr_i,wave=effec_wave[2]))/1.0e-17
    
    if photometry == 'Panstarr':
        effec_wave = [4810.,6170.,7520.,8660.]
        arr_f1[0],arr_f1_err[0] = np.array(mag2flux(psfMag_g,psfMagErr_g,wave=effec_wave[0]))/1.0e-17
        arr_f1[1],arr_f1_err[1] = np.array(mag2flux(psfMag_r,psfMagErr_r,wave=effec_wave[1]))/1.0e-17
        arr_f1[2],arr_f1_err[2] = np.array(mag2flux(psfMag_i,psfMagErr_i,wave=effec_wave[2]))/1.0e-17
    #arr_f1[3],arr_f1_err[3] = np.array(mag2flux(psfMag_z+0.040,psfMagErr_z,wave=effec_wave[3]))/1.0e-17
   	
    #curve_fit
    params, err = curve_fit(cmd,[SumFlux_g,SumFlux_r,SumFlux_i], arr_f1, sigma = arr_f1_err)

    index_blue = np.where(lam < 5900)
    index_red = np.where(lam > 5900)
    lam_red = lam[index_red]
    lam_blue = lam[index_blue]
    """
    flux1 = np.median(flux[wave_index1]*params[0])
    flux2 = np.median(flux[wave_index2]*params[1])
    flux_max = np.max([flux1,flux2])
    flux_min = np.min([flux1,flux2])
    """
    Res_lam = lam
    Res_flux = flux * params
    Res_ivar = ivar/(params[0]**2.)
    #Res_ivar = np.hstack((ivar[index_blue],ivar[index_red]))
    #Res_ivar = [ivar[index_blue]/(params[0]^2.),ivar[index_red]/(params[1]^2.)]
    index_ivar0 = np.where(Res_ivar == 0)
    #Res_ivar[index_ivar0] = 0.1
    index_ivar1 = np.where((Res_lam >5700)&(Res_lam < 5900))
    Res_ivar[index_ivar1] = 0.01
    Res_andmask = andmask
    Res_ormask  = ormask
    Res_num = len(Res_lam)
    Res_data    = np.zeros([5,Res_num])
    Res_data[0] = Res_flux
    Res_data[1] = Res_ivar
    Res_data[2] = Res_lam
    Res_data[3] = Res_andmask
    Res_data[4] = Res_ormask
    
    return Res_data, arr_f1, lam_blue, lam_red