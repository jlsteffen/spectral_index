# spectral_index

A module containing class and functions for performing astronomical operations. The package is designed to do the following;
- Identify objects from field images
- Extract fluxes from the identified objects
- Derive spectral indices

### Object Identification

The algorithm identifies object by searching for pixels above a given threshold in the image. The threshold in the function is given as a level of significance above the mean flux. A bitmask is created where pixels above the threshold are defined as 1 while the others are 0.

The algorithm then iteratively slides a window of a given size across the image. Within an iteration, if there is more than one pixel above the threshold, a centroiding function is applied to the frame to calculate the object center. The centroiding function is a simple center of mass type function which returns the x,y pixel coordinates of the object.

The coordinates are saved into a table that will be used for further operations. If a header with WCS information is given, the right ascension and declination of the objects are also given. Also a figure of the field with the identified objects marked will also be optionally produced. 

```
from astropy.io import fits
import numpy as np
from glob import glob
import pandas as pd
import spectral_index as spec
import os

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Open images
filepath = '.'

dat_files = glob(filepath + '/vlass_sample/*.tt0.subim.fits')

# Open minimum frequency image
im  = np.squeeze(fits.open(dat_files[7])[0].data)   
header = fits.open(dat_files[7])[0].header

# Run Object detection algorithm on the lowest Frequency image.

# Define a threshold to identify objects from.
# In this case I chose to set the threshold at 10 sigma over the image's mean value.
sigma = 10

# The algorithm will use a search window that iteratively passes over the image. 
# If pixels are found over the above threshold, a centroiding algorithm will be run to extract the object's center.
#
# The Search window should be large enough to capture an object, 
# but not so large that it will capture several objects in a single frame.
x_win, y_win = 200, 200

# The code will write a table with object positions to the given filepath.
# If a header with wcs information is provided, the code will calculate the RA and Dec of the objects.
# If fig = True, an image with the identified objects overlaid will be produced.

spec.findobj.finder(image = im, sigma = sigma, x_win = x_win, y_win = y_win, outpath = filepath, header=header, fig=True)
```

|	| objid	| xcoord	| ycoord |	ra	| dec |
|-|-------|---------|--------|------|-----|
|0|obj-1|	272.5590045|	1334.659359|	53.70636361|	-36.79323846|
|1|obj-2|	177.0944922|	5039.900904|	53.72140255|	-36.17560855|
|2|obj-3|	173.4299266|	-2473.983823|	53.73212263|	-37.42793586|
|3|obj-4|	727.7621614|	2378.348585|	53.6105068|	-36.61971755|
|4|obj-5|	880.2583103|	5499.500191|	53.57576913|	-36.09965396|



### Flux Extraction


### Spectral Indices

