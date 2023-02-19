# PyCali
- try to change

- Python version to re-calibrate the LAMOST spectra with the SDSS/PanSTARRS1 multi-band photometric data （see reference from [LAMOST Quasar Survey](https://ui.adsabs.harvard.edu/abs/2022arXiv221212876J/abstract)）


# pipeline

##  crossmatch with SDSS
- upload LAMOST QSO catalog

[sdss.china-vo casjob](http://sdss.china-vo.org/casjobs/MyDB.aspx)

[sdss casjob](http://skyserver.sdss.org/casjobs/MyDB.aspx)


- crossmatch and get Mag_g,Mag_r,Mag_i or [OIII]5007 

```
## 找附近的点
SELECT
  m.name AS name,m.radeg AS ra1, m.decdeg AS dec1, n.distance,
  o.ra AS ra2, o.dec AS dec2,o.petroMag_g AS Petro_g,o.petroMag_r AS Petro_r,o.petroMag_i AS Petro_i
  FROM MyDB.turnOff AS m
  OUTER APPLY dbo.fGetNearByObjAllEq( m.radeg, m.decdeg, 3.0/60.) AS n
  LEFT JOIN PhotoObj AS o ON n.objid=o.objid
```  


```
## 找最近的点
SELECT
  m.radeg AS ra1, m.decdeg AS dec1, n.distance,
  o.ra AS ra2, o.dec AS dec2,o.petroMag_g AS Petro_g,o.petroMag_r AS Petro_r,o.petroMag_i AS Petro_i
  FROM MyDB.DATA AS m
  OUTER APPLY dbo.fGetNearestObjEq(m.radeg, m.decdeg, 3.0/60.) AS n
  LEFT JOIN PhotoObj AS o ON n.objid=o.objid

```

- download data
```
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
```

## crossmatch with PanSTARRS1


# do calibration
- CaliFlux.py

- CaliFlux_with_SDSS_mag.py

- CaliFlux_with_Panstar_mag.py

- CaliFlux_with_SDSS_OIII.py






