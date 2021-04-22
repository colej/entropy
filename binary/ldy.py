"""
 Return limb-darkening and gravity-darkening coefficients

  The coefficients for a 4-parameter limb-darkening law and the gravity
 darkening coefficient, y, are interpolated from the tables provided by Claret
 and Bloemen (2011A&A...529A..75C) for ATLAS stellar atmosphere models and an
 assumed micro-turbulent velocity xi=2 km/s calculated using a least-squares
 fit. Also available are limb-darkening coefficients for the SDSS u', g', r',
 i' and z' bands from Claret (2004A&A...428.1001C), but no gravity darkening
 coefficient is available for these bands.


  To interpolate the limb-darkening and gravity-darkening coefficients at a
  given effective temperature (T_eff in K), surface gravity (log(g) in cgs
  units) and metallicity ([M/H] logarithmic metal abundance relative to
  solar), create an instance of the class LimbGravityDarkeningCoeffs for the
  desired photometric band and then call this instance with the parameters
  (T_eff, log(g), [M/H]). The first 4 values in the array returned are the
  limb-darkening coefficients a_1 .. a_4 and the final value is the
  gravity-darkening coefficient, y.  If no gravity-darkening coefficient is
  available then the fifth returned value will be None.

  If the input values of T_eff, log(g) or [M/H] are outside the tabulated
  range then all coefficients are returned as NaN.

   To see the names of the available photometric bands use the function
   ellc.ldy.list_bands(). For the SDSS bands, 'u_' = u', 'g_' = g', etc.

 Example
 -------

 >>> from ellc.ldy import LimbGravityDarkeningCoeffs, list_bands
 ['B', 'C', 'H', 'I', 'J', 'K', 'Kp', 'R', 'S1', 'S2', 'S3', 'S4', 'U', 'V',
 >>> print(list_bands())
 'b', 'u', 'v', 'y', 'u_', 'g_', 'r_', 'i_', 'z_']
 >>> ldy_Kp = LimbGravityDarkeningCoeffs('Kp')
 >>> Teff = 10450
 >>> logg = 3.9
 >>> M_H = -0.31
 >>> a1,a2,a3,a4,y = ldy_Kp(Teff, logg, M_H)

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.io import fits
from os.path import join,abspath,dirname
from scipy.interpolate import RegularGridInterpolator

dir_path = '/lhome/colej/astro/ebs/coefficients'
ldc_filebase = 'g_claret2017_{}_LDC_{}.fits'
gdc_filebase = 'g_claret2017_{}_GDC.fits'
ldcKp_filebase = 'g_claret2011_{}_LDC_{}.fits'
gdcKp_filebase = 'g_claret2011_{}_GDC.fits'
combined_filebase = 'g_claret2017_{}_LDC_{}_GDC.fits'
ldc2011_filebase = 'g_claret2011_LDC_{}.fits'
gdc2011_filebase = 'g_claret2011_GDC.fits'

class LimbGravityDarkeningCoeffs:

    def __init__(self, band, law):
        self._band = band
        self._law = law
        print('Band: {}'.format(self._band))
        print('Law: {}'.format(self._law))

        # if band in ['Kp','Tp']:
        if band in ['Tp']:
        # if ((self._band == 'Tp') or (self._band=='Kp')) :

            # LDC tables are linear Teff
            ldc_hdulist = fits.open(join(dir_path,
                                ldc_filebase.format(self._band,self._law)))
            # GDC tables are log Teff
            gdc_hdulist = fits.open(join(dir_path,
                                gdc_filebase.format(self._band)))

            ldc_grid = ldc_hdulist[0].data
            gdc_grid = gdc_hdulist[0].data
            print(gdc_grid)

        elif band in ['Kp']:
        # if ((self._band == 'Tp') or (self._band=='Kp')) :

            # LDC tables are linear Teff
            ldc_hdulist = fits.open(join(dir_path,
                                ldcKp_filebase.format(self._band,self._law)))
            # GDC tables are log Teff
            gdc_hdulist = fits.open(join(dir_path,
                                gdcKp_filebase.format(self._band)))

            ldc_grid = ldc_hdulist[0].data#)[:,0,:,:,:]
            print(ldc_grid[1,:,:,:])
            gdc_grid = gdc_hdulist[0].data

        else:
            # LDC tables are linear Teff
            ldc_hdulist = fits.open(join(dir_path, ldc2011_filebase.format(self._law)))
            # GDC tables are log Teff
            gdc_hdulist = fits.open(join(dir_path, gdc2011_filebase))

            ldc_bands = ldc_hdulist['Band'].data
            gdc_bands = gdc_hdulist['Band'].data
            i_ldc_band = (ldc_bands['Band'].rfind(self._band) > -1).nonzero()[0][0]
            i_gdc_band = (gdc_bands['Band'].rfind(self._band) > -1).nonzero()[0][0]
            # print(i_ldc_band,i_gdc_band)
            # print(np.shape(ldc_hdulist[0].data))
            ldc_grid = (ldc_hdulist[0].data)[i_ldc_band,:,:,:,:]
            gdc_grid = (gdc_hdulist[0].data)[i_gdc_band,:,:,:,:]


        logTeff = gdc_hdulist['logTeff'].data
        teff = ldc_hdulist['Teff'].data

        logg_ldc = ldc_hdulist['logg'].data
        M_H_ldc = ldc_hdulist['M_H'].data

        logg_gdc = gdc_hdulist['logg'].data
        M_H_gdc = gdc_hdulist['M_H'].data

        ldc_hdulist.close()
        gdc_hdulist.close()

        ldc_pts = (M_H_ldc,logg_ldc,teff)
        gdc_pts = (M_H_gdc,logg_gdc,logTeff)

        self._y = RegularGridInterpolator(gdc_pts, gdc_grid, fill_value=None,method='linear',bounds_error=True)
        self._a1 = RegularGridInterpolator(ldc_pts,ldc_grid[0,:,:,:],fill_value=None)
        self._a2 = RegularGridInterpolator(ldc_pts,ldc_grid[1,:,:,:],fill_value=None)



    def __call__(self, teff, logg, M_H):
        if teff < 0:
            y,a1,a2  = [np.nan,np.nan,np.nan]
        else:
            try:
                y  = self._y ([M_H,logg,np.log10(teff)])
                a1 = self._a1([M_H,logg,teff])
                a2 = self._a2([M_H,logg,teff])
            except:
                y,a1,a2  = [np.nan,np.nan,np.nan]



        return (np.array([a1, a2, y])).flatten()


class LimbGravityDarkeningCoeffs_quadlaw:

    def __init__(self, band, law):
        self._band = band
        self._law = law


        coeff_hdulist = fits.open(join(dir_path,
                            combined_filebase.format(self._band,self._law)))
        teff = coeff_hdulist['Teff'].data

        logg = coeff_hdulist['logg'].data
        M_H = coeff_hdulist['M_H'].data
        # Xi = ldc_hdulist['Xi'].data

        coeff_grid = coeff_hdulist[0].data

        coeff_hdulist.close()

        coeff_pts = (M_H,logg,teff)

        self._a1 = RegularGridInterpolator(coeff_pts,coeff_grid[0,:,:,:],
                                               fill_value=None)
        self._a2 = RegularGridInterpolator(coeff_pts,coeff_grid[1,:,:,:],
                                               fill_value=None)
        self._a3 = RegularGridInterpolator(coeff_pts,coeff_grid[2,:,:,:],
                                               fill_value=None)
        self._a4 = RegularGridInterpolator(coeff_pts,coeff_grid[3,:,:,:],
                                               fill_value=None)
        self._y = RegularGridInterpolator(coeff_pts,coeff_grid[4,:,:,:],
                                               fill_value=None)



    def __call__(self, teff, logg, M_H):
        if teff < 0:
            y,a1,a2  = [np.nan,np.nan,np.nan]
        else:
            try:
                y  = self._y ([M_H,logg,teff])
                a1 = self._a1([M_H,logg,teff])
                a2 = self._a2([M_H,logg,teff])
                a3 = self._a3([M_H,logg,teff])
                a4 = self._a4([M_H,logg,teff])
            except:
                y,a1,a2,a3,a4  = [np.nan,np.nan,np.nan,np.nan,np.nan]



        return (np.array([a1, a2, a3, a4, y])).flatten()
