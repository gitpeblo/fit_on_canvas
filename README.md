fit_on_canvas.py
================

/ DESCRIPTION /-----------------------------------------------------------------

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

Instead of the IMAGE (x,y) pixel positions used to center each input image, the
user can provide a single WCS (RA,Dec) coordinate, which will be converted into
pixel coordinates using the header information (useful if the images are
astrometrically aligned).

A buffer (square annulus) can be placed around the images when they are placed
on their canvases.
If BUFFER is set to "None", the FOVs of the images will be projected on the
canvases as they are (and no BUFFER is considered, i.e. it is set to 0).
Otherwise, the images are first trimmed to minimize their size (i.e. in order to
only keep the valid data points), and then the desired BUFFER is added.

/ USAGE / ----------------------------------------------------------------------

    fit_on_canvas.py -h

PRODUCT: canvas_\<fits\>, canvas_\<fits2\>, .. <BR>
RETURNS: list of output files     

EXAMPLES:

    fit_on_canvas.py -FILE image1.fits -FILE image2.fits -WCS 62.466021 -56.118515 -BUFFER 0

   Will create "canvas\_image1.fits" and "canvas\_image2.fits", with each image
   centered around (RA, Dec) = (62.466021, -56.118515).
   The images are trimmed to have no borders (BUFFER = 0), therefore at least
   one of the output images will "touch" the canvas border.


