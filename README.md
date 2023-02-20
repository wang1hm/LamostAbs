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



# download LAMOST spectra

[LAMOST spectra](http://www.lamost.org/lmusers/)
LAMOST的光谱数据在官网上下载的。
选择某个DR，然后选择low resolution

[LAMOST dr9](http://www.lamost.org/dr9/v1.1/search)
这里就可以直接用ra，dec进行下载了。


# download SDSS spectra
[dr16.sdss](https://dr16.sdss.org/optical/spectrum/search)

[dr17.sdss](https://dr17.sdss.org/optical/spectrum/search)


# do calibration in prep

LAMOST光谱流量定标

== 使用测光数据进行流量定标：1）由于LAMOST的波长覆盖范围有限，所以只用到了gri三个波段，这三个波段可以在SDSS或者Pan-starr上获得。2）用SDSS（或者Pan-starr）的透过率曲线去卷积LAMOST光谱，就得到了对应的单色光度，用这个光谱去拟合在SDSS（或者Pan-starr）上得到的实际光度，就可以实现流量定标。

== 利用窄发射线做流量定标:这种方法是假设AGN的某些窄发射线是不发生变化的，将某个光谱乘上一个值，是的其窄线的flux和另一条光谱一致，就完成了流量定标，一般用到的窄线是[OIII]5007.

- CaliFlux.py

- CaliFlux_with_SDSS_mag.py

- CaliFlux_with_PanSTAR_mag.py

- CaliFlux_with_SDSS_OIII.py






