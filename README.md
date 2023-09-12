# spectral_index

A module containing class and functions for performing astronomical operations. The package is designed to do the following;
- Identify objects from field images
- Extract fluxes from the identified objects
- Derive spectral indices

### Object Identification

The algorithm identifies object by searching for pixels above a given threshold in the image. The threshold in the function is given as a level of significance above the mean flux. A bitmask is created where pixels above the threshold are defined as 1 while the others are 0.

The algorithm then iteratively slides a window of a given size across the image. Within an iteration, if there is more than one pixel above the threshold, a centroiding function is applied to the frame to calculate the object center. The centroiding function is a simple center of mass type function which returns the x,y pixel coordinates of the object.

The coordinates are saved into a table that will be used for further operations. If a header with WCS information is given, the right ascension and declination of the objects are also given. Also a figure of the field with the identified objects marked will also be optionally produced. 

|	| objid	| xcoord	| ycoord |	ra	| dec |
|-|-------|---------|--------|------|-----|
|0|obj-1|	272.5590045|	1334.659359|	53.70636361|	-36.79323846|
|1|obj-2|	177.0944922|	5039.900904|	53.72140255|	-36.17560855|
|2|obj-3|	173.4299266|	-2473.983823|	53.73212263|	-37.42793586|
|3|obj-4|	727.7621614|	2378.348585|	53.6105068|	-36.61971755|
|4|obj-5|	880.2583103|	5499.500191|	53.57576913|	-36.09965396|



### Flux Extraction


### Spectral Indices

