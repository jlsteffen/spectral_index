# spectral_index

A module containing class and functions for performing astronomical operations. The package is designed to do the following;
- Identify objects from field images
- Extract fluxes from the identified objects
- Derive spectral indices

### Object Identification

The algorithm identifies object by searching for pixels above a given threshold in the image. The threshold in the function is given as a level of significance above the mean flux. A bitmask is created where pixels above the threshold are defined as 1 while the others are 0.

The algorithm then iteratively slides a window of a given size across the image. Within an iteration, if there is more than one pixel above the threshold, a centroiding function is applied to the frame to calculate the object center. The centroiding function is a simple center of mass type function which returns the x,y pixel coordinates of the object.

The coordinates are saved into a table that will be used for further operations. If a header with WCS information is given, the right ascension and declination of the objects are also given. Also a figure of the field with the identified objects marked will also be optionally produced. 

An example of the object identification is shown below;
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

The script outputs the following table from the vlass sample.

|	| objid	| xcoord	| ycoord |	ra	| dec |
|-|-------|---------|--------|------|-----|
|0|obj-1|	272.5590045|	1334.659359|	53.70636361|	-36.79323846|
|1|obj-2|	177.0944922|	5039.900904|	53.72140255|	-36.17560855|
|2|obj-3|	173.4299266|	-2473.983823|	53.73212263|	-37.42793586|
|3|obj-4|	727.7621614|	2378.348585|	53.6105068|	-36.61971755|
|4|obj-5|	880.2583103|	5499.500191|	53.57576913|	-36.09965396|
|.|...|...|...|...|...|

If the figure function is enabled, the following figure is produced highlighting the object identification results.

![Figure 1](https://github.com/jlsteffen/spectral_index/blob/main/figs/field.png)

### Flux Extraction

After a table of object is identified, the module can be used to extract fluxes for the objects. The function uses a circular aperture of a given radius. The function loops over each frequency frame given and each object identified and extracts the object's flux. The flux and frequency for each object is then saved into a table so that the spectral index can later by calculated.

An example of the flux extraction is shown below;
```
# Extract fluxes from the identified objects.

# The code will loop over each frequency frame and each object identified and extract the flux
# for each object at each frequency. 

# first define a list of files to loop over. Each file should be an image taken at a different frequency.
files = dat_files

# The provide the table of the objects with their coordinates
table = pd.read_csv('objects.csv', delimiter=',', index_col=0)

# Define the radius of the circular aperture to extract the fluxes from
radius = 50

# The code will write tables containing the frequencies and fluxes for each object.
# Provide a filepath to save these table to.
outpath = filepath + '/objects'
if not os.path.exists(outpath):
    os.makedirs(outpath)

spec.extract_flux.extract_objects(files = files, table = table, radius = radius, outpath = outpath)
```
The script produces a table like the one below for each identified object.

| |      frequency |         flux |
|-|----------------|--------------|
|0|   3.308007e+09 |  6451.352935 |
|1|   2.412005e+09 | 17333.243388 |
|2|   2.156004e+09 | 23178.419305 |
|3|  2.796006e+09 | 12489.795790 |
|4|   3.564007e+09 |  5578.265733 |
|-|----------------|--------------|

### Spectral Indices

