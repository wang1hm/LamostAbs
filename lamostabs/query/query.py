#!/usr/bin/env python
from astropy import units as u
from astropy.table import join
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
import os, warnings
import pandas as pd


class SDSSQuery():
    """_summary_
    """

    def __init__(self, service='SDSS') -> None:
        """_summary_

        Parameters
        ----------
        service : str, optional
            _description_, by default 'SDSS'

        Raises
        ------
        NotImplementedError
            _description_
        """
        service = service.upper()
        if 'CDS' in service:
            service = 'CDS_XMATCH'
        try:
            if service not in ['SDSS', 'CDS_XMATCH']:
                raise NotImplementedError("The input service is not in implemented! Select one from ['SDSS', 'CDS_XMATCH'].")
        except NotImplementedError as exp:
          warnings.warn(exp, "Fall back to `service = 'SDSS'`.")
        self.service = service


    def upload_match(self, table_in, ra_col='ra', dec_col='dec', radius=2.0 * u.arcsec):
        """_summary_

        Parameters
        ----------
        table_in : astropy.table.Table object
            _description_
        ra_col : str, optional
            _description_, by default 'ra'
        dec_col : str, optional
            _description_, by default 'dec'
        radius : astropy Quantity, optional
            _description_, by default 2.0*u.arcsec

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        KeyError
            _description_
        """
        tb = Table(table_in)
        old_cols = tb.columns
        new_cols = [x.lower() for x in old_cols]
        tb.rename_columns(old_cols, new_cols)
        try:
            ra = tb[ra_col]
            dec = tb[dec_col]
        except KeyError as e:
            raise KeyError("Column name(s) for ra and/or dec not found: {}, {}.".format(ra_col, dec_col))
        try:
            coord = SkyCoord(ra=ra.value*u.deg, dec=dec.value*u.deg, frame='icrs')
        except TypeError as e_t:
            ra_s = pd.Series(ra.value)
            dec_s = pd.Series(dec.value)
            co_str = ra_s + ' ' + dec_s
            coord = SkyCoord(co_str, unit=(u.hourangle, u.deg), frame='icrs')
        if self.service == 'SDSS':
            from astroquery.sdss import SDSS
            res = SDSS.query_crossid_async(coordinates=coord, radius=radius)
            res_tb = Table.read(res.text, format='ascii')
            table_out = join(left=tb, right=res_tb)
        elif self.service == 'CDS_XMATCH':
            table = XMatch.query_async(cat1=tb,
                                       cat2='vizier:V/147/sdss12',
                                       max_distance=radius, 
                                       colRA1=ra_col,
                                       colDec1=dec_col)
            table_out = table[table['mode']==1]
        self.table_out = table_out
        return table_out
        