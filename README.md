fit_on_canvas
=============

Given an input set of FITS images and corresponding (x,y) pixel positions
("centers"), this program expands each input image in order to accommodate all of
them within a larger, common canvas.
The field-of-views (FOVs) of the input images are expanded such as to keep the
input (x,y) pixel positions at the center of the new canvas.

First, the common canvas size is determined adopting the maximum distance between
the (x,y) centers and the border of each mage.
Then, each image is re-plotted on the larger canvas (the number of output images
corresponds to the number of input images).

The WCS coordinates of the output image are kept in place by updating the CRPIX
keywords.

A buffer (square annulus) can be placed around the images when they are placed
on their canvases.
If BUFFER is set to "None", the images will be projected on the canvas as they
are (and no BUFFER is considered, i.e. it is set to 0).
Otherwise, the images are first trimmed to minimize their size (i.e. in order to
only keep the valid data points), and then the desired BUFFER is added.
