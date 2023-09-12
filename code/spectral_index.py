#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is designed to detect sources within astronomical images, extract
fluxes from the identified sources, and derive spectral indices for the 
identified objects.
"""

import numpy as np
import math as m
import pandas as pd
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from copy import copy


def nanlog10(x):
    """
    Calculates the log10 of an array and handles for invalid values. Invalid
    values ( x <= 0 ) are calculated as NaNs.

    Parameters
    ----------
    x : n-d array
        Array to calculate log10 on.

    Returns
    -------
    values : n-d array
        Log10 of the input array.

    """
    values = np.log10(x, out=np.zeros(x.shape) * np.nan, where=x>0)
    return values

#_____________________________________________________________________________

class findobj:
    """
    Class to run the object identification and centroiding algorithm.
    """
    def __init__(self):
        self.centroid()
        self.finder()
        self.field_img()
        
    def finder(image, sigma, x_win, y_win, outpath, header=False, fig=False):
        """
        Detects sources within an image based on a given sigma level. The code 
        defines a threshold on which to identify objects based on a multiple of 
        the standard deviation. Pixels below the threshold are set to 0 and pixels
        above the threshold are set to 1. 
        
        The function then takes a defined rectangular search window which 
        iteratively loops over the image. Windows with more than 1 pixel over the
        given threshold will calculate a centroid. Empty windows are passed over.
        
        An array of centroids is returned.
        
    
        Parameters
        ----------
        image : 2D-array
            Science Image to detect sources on.
        sigma : float
            The sigma level to detect objects on (e.g. sigma=5 will detect sources
            5-sigma over the image's mean).
        x_win : float
            x dimension size of the search window.
        y_win : float
            y dimension size of the search window.
        outpath : str
            Filepath to save the object table and figure if requested.
        header : fits header, optional
            Optional. Header of image. Will use header to extract WCS 
            information which is used to calculate the right ascension and 
            declination. If false, no sky coordinates will be added to the 
            object file. The default is False.
        fig : bool, optional
            Toggle to plot a field showing the object detection. The default 
            is False.
    
        Returns
        -------
        None.
    
        """
        
        threshold = np.mean(image) + (sigma * np.std(image))
        
        # Creates a bit image where 0 = empty, 1 = source
        bitim = 0 * (image <= threshold) + 1 * (image > threshold)
        
        x_cents = []
        y_cents = []
        for i in range( m.ceil(image.shape[1] / x_win) ):
            for j in range( m.ceil(image.shape[0] / y_win) ):
                bitwin = bitim[j*y_win : (j+1)*y_win , i*x_win : (i+1)*x_win]
                
                # skip frames with less than or equal to one pixel above the 
                # threshold
                if np.nansum(bitwin) <= 1:
                    pass
                else:
                    # Find centroid
                    x_cent, y_cent = findobj.centroid(bitwin)
                    
                    # Centroid is just of the given sub_frame
                    # Need to get the image coordinates
                    x_cents.append(x_cent + i * x_win)
                    y_cents.append(y_cent + j * y_win)
        x_cents = np.array(x_cents)
        y_cents = np.array(y_cents)

        # Generate obj names
        objid = ['obj-'+str(i+1) for i in range(len(x_cents))]

        # Make table        
        DF = pd.DataFrame(np.array([objid,x_cents,y_cents]).T, columns=['objid','xcoord','ycoord'], index=None)
        
        # Read pix/world from header
        if header != False:
            # extract wcs from header
            w = WCS(header, naxis=2)
            
            sky = w.pixel_to_world(x_cents,y_cents)
            
            ra = sky.ra.value
            dec = sky.dec.value
            DF['ra'] = ra
            DF['dec'] = dec
            
        # Write csv
        DF.to_csv(outpath+'/objects.csv')
        
        # print figure if asked
        findobj.field_img(image=image, xcoord=x_cents, ycoord=y_cents, radius=50, outpath = outpath+'/field.png')
        
    def centroid(data):
        """
        Finds the centroid of a given rectangular image/image cutout.
        
        Parameters
        ----------
        data : 2d-array
            2d image array
        
        Returns
        -------
        x_cent : float
            x centroid coordinate
        
        y_cent : float
            y centroid coordinate
        """
        
        x_grid = np.arange(0, data.shape[1], 1)
        y_grid = np.arange(0, data.shape[0], 1)
        
        x_cent = np.sum( np.dot(data, x_grid) ) / np.sum( data )
        y_cent = np.sum( np.dot(y_grid, data) ) / np.sum( data )
        
        return x_cent, y_cent

    
    def field_img(image, xcoord, ycoord, radius, outpath):
        """
        Prints the field image with the identified objects labeled with the 
        apertures from which the fluxes will be extracted.
    
        Parameters
        ----------
        image : 2D-array
            Science Image to detect sources on.
        xcoord : 1-d array
            x-coords of the object center.
        ycoord : 1-d array
            y-coords of the object center..
        radius : float
            Radius of the extraction apertures.
        outpath : str
            filepath to save the figure.
    
        Returns
        -------
        None.
    
        """
        f, ax = plt.subplots()
        
        vmin, vmax = 0, image.max() / 10
        plt.imshow(image, cmap='Greys_r', vmin=vmin, vmax=vmax)
        
        for i in range(len(xcoord)):
            psf_circle = patches.Circle([xcoord[i], ycoord[i]], radius=radius, ec='r', fc='none')
            ax.add_patch(copy(psf_circle))
        
        ax.xaxis.set_tick_params(which ='both', top=True, direction='in', labelbottom=False, labeltop=False)
        ax.yaxis.set_tick_params(which ='both', right=True, direction='in', labelleft=False, labelright=False)
        
        plt.savefig(outpath, bbox_inches='tight')
        plt.close(f)

#_____________________________________________________________________________
class extract_flux:
    """
    Class to extract fluxes from the identified objects.
    """
    def __init__(self):
        self.extract_objects()
        self.aperture()
        self.flux_extract()

    def extract_objects(files, table, radius, outpath):
        """
        Iterates over each identified object and each frequency band to generate
        a table of fluxes and frequencies. Writes a csv table for each object 
        containing fluxes and frequencies.
    
        Parameters
        ----------
        files : list
            list of filepaths for the files to be operated on.
        table : pandas dataframe
            dataframe containing object positions of the identified sources. The
            table should contain objid, xcoord, and ycoord columns.
        radius : float
            radius of the aperture that the fluxes are extracted from.
        outpath : str
            file path to which the object tables are save to.
    
        Returns
        -------
        None.
    
        """        
        # load in all frames first for efficiency
        frames = []
        freq = []
        for i in range(len(files)):
            fil = fits.open(files[i])[0]
            frames.append( np.squeeze(fil.data) )
            freq.append( fil.header['restfrq'] )    
        
        
        for i in range(len(table)):    
            flux = []
            for j in range(len(frames)):
                im  = frames[j]
                
                #t1 = t.time()
                aper = extract_flux.aperture(im, table['xcoord'][i], table['ycoord'][i], radius)
                #print('aperture time = '+str(t.time() - t1))
                
                flux.append( extract_flux.flux_extract(im, aper) )
                
            flux = np.array(flux)
            
            DF = pd.DataFrame(np.array([freq, flux]).T, columns=['frequency', 'flux'], index=None)
            
            DF.to_csv(outpath+'/'+table['objid'][i]+'.csv')
            print('\nTable written to '+outpath+'/'+table['objid'][i]+'.csv'+'\n')  
    
    def aperture(im, x, y, radius):
        """
        Creates a bitmask for a circular aperture selection. 
    
        Parameters
        ----------
        im : 2-D array
            The image/array from which to apply the aperture on.
        x : Float
            x-coord of the object center.
        y : Float
            y-coord of the object center.
        radius : Float
            Radius of the aperture.
    
        Returns
        -------
        mask : bitmask defining a circular aperture.
    
        """
        # setup meshgrid
        y_grid,x_grid = np.ogrid[0: im.shape[0], 0: im.shape[1]]
        
        # center grid on object
        x_cent = x_grid - x
        y_cent = y_grid - y
        
        # create mask: TRUE = pixels within the aperture
        mask = x_cent**2 + y_cent**2 <= radius**2
        
        return mask
    
    def flux_extract(im, mask):
        """
        Extracts flux from a given image through a given aperture.
    
        Parameters
        ----------
        im : 2-D array
            The image/array from which to apply the aperture on.
        mask : 2-D array
            mask array giving the aperture to select flux from.
    
        Returns
        -------
        flux : Flux density.
    
        """
        extract = im * mask
        
        flux = np.nansum(extract) * np.nansum(mask)
        
        return flux

#_____________________________________________________________________________
class derive_specind:
    """
    Class to derive spectral indices with the extracted fluxes and given
    frequencies.
    """
    def __init__(self):
        self.derive()
        self.specind()
        self.spectral_fig()

    def derive(tablepath, objectpath, fig=False):
        """
        Add a spectral index column to the object table. Calculates the spectral
        index for all objects in the table based on their frequency and extracted
        fluxes.
    
        Parameters
        ----------
        filename : str
            Filepath of the Master table containing the object ID's and 
            positions.
        objectpath : str
                Filepath to the individual object tables.
        fig : bool, optional
            Toggle to plot a SED for each object. The default is False.
    
        Returns
        -------
        None.
    
        """
        table = pd.read_csv(tablepath, delimiter=',', index_col=0)
        
        spectral = []
        for i in range(len(table)):
            ff = pd.read_csv(objectpath+table['objid'][i]+'.csv', delimiter=',', index_col=0)
            
            alpha, b = derive_specind.specind(ff['frequency'], ff['flux'])
            
            spectral.append(alpha)
            
            if fig==True:
                derive_specind.spectral_fig(ff['frequency'], ff['flux'], 
                        table['objid'][i], objectpath+table['objid'][i]+'.png')
            
        spectral = np.array(spectral)
        
        table['alpha'] = spectral
        
        table.to_csv(tablepath)
        
        print('\nSpectral indices added to '+tablepath)
    
    def specind(frequency, flux):
        """
        Solves the spectral index for a source given a set of frequencies and 
        fluxes. The spectral index is solved by taking the slope between the log 
        of the flux and the log of the frequency.
        
        log ( flux ) = alpha * log ( frequency )
    
        Parameters
        ----------
        frequency : 1-d array
            An array of frequencies.
        flux : 1-d array
            An array of fluxes.
    
        Returns
        -------
        alpha : Float
            The spectral index.
        b : Float
            A constant, the y-intercept of the linear equation
    
        """
        
        logfreq = nanlog10(frequency)
        logflux = nanlog10(flux)
        
        # mask NaN values
        mask = (np.isnan(logfreq))|(np.isnan(logflux))
        
        alpha, b = np.polyfit(logfreq[~mask], logflux[~mask], 1)
        
        return alpha, b
    
    def spectral_fig(frequency, flux, objid, outpath):
        """
        Prints figures of the object's SED with spectral index derivation.
    
        Parameters
        ----------
        frequency : 1-d array
            An array of frequencies.
        flux : 1-d array
            An array of fluxes.
        objid : STR
            objid or the name of the object.
    
        Returns
        -------
        None.
    
        """
    
        f, ax = plt.subplots()
        
        plt.scatter(nanlog10(frequency), nanlog10(flux), c='r', marker='D')
        
        alpha, b = derive_specind.specind(frequency, flux)
        
        xx = np.linspace(9.3,9.6, 100)
        yy = alpha * xx + b
        plt.plot(xx, yy, c='k', linestyle='--')
        
        fontsize = 16
        plt.title(objid, fontsize=fontsize)
        plt.xlabel(r'log(Frequency, Hz)', fontsize=fontsize)
        plt.ylabel(r'log(Flux Density, Jy)', fontsize=fontsize)
        plt.text(.99, .99, r'$\alpha$ $=$ '+str(alpha), ha='right', va='top', 
                 transform=ax.transAxes, fontsize=fontsize/1.25)
        f.tight_layout()
        
        plt.savefig(outpath, bbox_inches='tight')
        plt.close(f)