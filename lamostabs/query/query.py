#!/usr/bin/env python
from astropy import units as u
from astropy.table import join
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
import os, warnings
import pandas as pd
import numpy as np


class SDSSQuery():
    """Class for querying the SDSS photometric catalog.
    """

    def __init__(self, service='SDSS') -> None:
        """Initialize an instance of SDSSQuery.

        Parameters
        ----------
        service : str, optional
            Web API used for the query, can be 'SDSS' (based on `astroquery.sdss.SDSS`) or 'CDS_XMATCH' (based on `astroquery.xmatch.XMatch`), by default 'SDSS'

        Raises
        ------
        NotImplementedError
            Currently only 'SDSS' and 'CDS_XMATCH' are implemented.
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
        """Upload and crossmatch the local table with SDSS photometric catalog using coordinates.

        Parameters
        ----------
        table_in : astropy.table.Table object
            The input table to be uploaded
        ra_col : str, optional
            The column name of R.A. in `table_in`, by default 'ra'
        dec_col : str, optional
            The column name of Decl. in `table_in`, by default 'dec'
        radius : astropy Quantity, optional
            Radius for the crossmatch, by default 2.0*u.arcsec

        Returns
        -------
        astropy.table.Table object
            The joined table after crossmatch

        Raises
        ------
        KeyError
            Raised when `ra_col`/`dec_col` is not found in the column names of the input table.
        """
        tb = Table(table_in)
        tb['upload_id'] = np.arange(len(tb))
        old_cols = tb.colnames
        new_cols = [x.lower() for x in old_cols]
        tb.rename_columns(old_cols, new_cols)
        try:
            ra = tb[ra_col]
            dec = tb[dec_col]
        except KeyError as e:
            raise KeyError("Column name(s) for ra and/or dec not found: {}, {}.".format(ra_col, dec_col))
        try:
            coord = SkyCoord(ra=ra.data*u.deg, dec=dec.data*u.deg, frame='icrs')
        except TypeError as e_t:
            ra_s = pd.Series(ra.data)
            dec_s = pd.Series(dec.data)
            co_str = ra_s + ' ' + dec_s
            coord = SkyCoord(co_str, unit=(u.hourangle, u.deg), frame='icrs')
        if self.service == 'SDSS':
            from astroquery.sdss import SDSS
            res = SDSS.query_crossid_async(coordinates=coord, radius=radius, obj_names=tb['upload_id'])
            res_tb = Table.read(res.text, format='ascii')
            table_out = join(left=tb, right=res_tb, keys_left='upload_id', keys_right='obj_id')
        elif self.service == 'CDS_XMATCH':
            res = XMatch.query_async(cat1=tb,
                                       cat2='vizier:V/147/sdss12',
                                       max_distance=radius, 
                                       colRA1=ra_col,
                                       colDec1=dec_col)
            table = Table.read(res.text, format='ascii')
            table_out = table[table['mode']==1]
        self.table_out = table_out
        return table_out
        