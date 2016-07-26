#!/usr/bin/python

'''
/ DESCRIPTION / ----------------------------------------------------------------

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

USAGE: fit_on_canvas.py -h
PRODUCT: canvas_<fits1>, canvas_<fits2>, ..
RETURNS: list of output files     

EXAMPLES:

$> fit_on_canvas.py -FILE image1.fits -FILE image2.fits -WCS 62.466021 -56.118515 -BUFFER 0

   Will create "canvas_image1.fits" and "canvas_image2.fits", with each image
   centered around (RA, Dec) = (62.466021, -56.118515).
   The images are trimmed to have no borders (BUFFER = 0), therefore at least
   one of the output images will "touch" the canvas border.

/ HISTORY / --------------------------------------------------------------------

26/07/2016: fit_on_canvas_v0.py /

> Update:
> Corrected:
> Fixed:

/ MEMO / -----------------------------------------------------------------------

/ NOTICE / ---------------------------------------------------------------------

--------------------------------------------------------------------------------
'''

import pyfits
import sys
import os
import argparse
import re
import numpy as np
import warnings
import pywcs

################################################################################
def fit_on_canvas(fnames, centers, coord_type, **optional_parameters):
  """
  Adapts input images onto a larger canvas that accommodates all input images.
  The center of the input image in the output canvas is specified by the user.
  
  USAGE: fit_on_canvas(FITS_names[list],centers_x_y[array,list_of_array], coord_type["WCS","IMAGE"][, BUFFER[pixel,None], SILENT[True,False]])
  PRODUCT: canvas_<fits1>, canvas_<fits2>, ..
  RETURNS: list_of_canvas  
  """

  #-----------------------------------------------------------------------------
  ADD_BUFFER = False
  BUFFER     = None

  # If buffer is set to "None", the images will be projected on the canvas as
  #   they are (and no BUFFER is considered, i.e. it is set to 0)
  # Otherwise, the images are first trimmed to minimize their size (i.e. in
  #   order to only keep the valid data points), and then the desired BUFFER is
  #   added
  
  if ('BUFFER' in optional_parameters):
    if optional_parameters['BUFFER'] is not None:
      APPLY_BUFFER = True
      BUFFER       = int(optional_parameters['BUFFER'])
    else:
      APPLY_BUFFER = False
      BUFFER       = 0
      
    
  SILENT = False
  
  if ('SILENT' in optional_parameters):
    if (optional_parameters['SILENT'] != False):
      SILENT = optional_parameters['SILENT']
  #-----------------------------------------------------------------------------

  coords = []
  # array of actual (x,y) centering coordinates
  # (i.e., input IMAGE or converted from input WCS to IMAGE)
  
  if ( (coord_type == "IMAGE" ) & (len(fnames) != len(centers)) ):
    print "ERROR: Number of input files doesn't match number of input centers!"
    exit

  fits_names = []
  # file names (stripped of extension)
  fits_exts = [0]*len(fnames)
  # FITS extensions of input files [default=0]

  for f, fname in enumerate(fnames):
    match = re.search('(.*fits)(\[(\d+)\])?', fname)
    fits_name = match.group(1)
    fits_ext  = match.group(3)

    fits_names.append(fits_name)
    fits_exts.append(fits_ext)


  if (coord_type == "IMAGE" ):
  
    coords = [( int(round(centers[f][0])) , int(round(centers[f][1])) ) for f in range(len(centers))]
 
  if (coord_type == "WCS"):

    for f, fits in enumerate(fits_names):
      HDU_fits = pyfits.open(fits)
      image  = HDU_fits[fits_exts[f]].data
      header = HDU_fits[fits_exts[f]].header
      wcs = pywcs.WCS(header)

      # Converting WCS coordinates to IMAGE coordinates:
      sky_coord = np.array(centers, np.float_)
      # converting centers (input WCS) to numpy array for wcs.wcs_sky2pix
      pix_coords = wcs.wcs_sky2pix(sky_coord, 1)
      # converting WCS coordinates to IMAGE coordinate
      
      coords.append( pix_coords[-1].tolist() )
      coords = [( int(round(coords[f][0])) , int(round(coords[f][1])) ) for f in range(len(coords))]


  if not (SILENT):
    print "Input parameters:"

    for f, fits in enumerate(fits_names):
      print "  %s [%s] -> centering with BUFFER = %s @ image (%s,%s)" % (fits_names[f],fits_exts[f],BUFFER,coords[f][0],coords[f][1])

  #-----------------------------------------------------------------------------

  
  # defining canvas size -------------------------------------------------------
  fits_NAXIS1 = []
  fits_NAXIS2 = []
  # arrays containing NAXIS1 and NAXIS2 keywords of all files
  
  # min/MAX x and y coordinates of image:
  
  # Could be either of NAXIS* or the value of the "extremes" of the image
  # The "extremes" of the image are defined as the left-most, right-most, top,
  #   or bottom limits of the actual data image (i.e. excluding the pixels
  #   which are just filling the canvas, and which usually are set to 0)
  
  images_x_min = []
  images_y_min = []
  images_x_MAX = []
  images_y_MAX = []
  # [pixel]
  

  for f, fits in enumerate(fits_names):
    HDU_fits = pyfits.open(fits)
    image  = HDU_fits[fits_exts[f]].data
    header = HDU_fits[fits_exts[f]].header
    fits_NAXIS1.append( HDU_fits[fits_exts[f]].header['NAXIS1'] )
    fits_NAXIS2.append( HDU_fits[fits_exts[f]].header['NAXIS2'] )

    # Coordinates of extreme pixels:
    #
    # Here, np.where returns the <tuples> with the indexes corresponding to
    #   valid (i.e. > 0) pixels
    # The first array (<index> = 0) created by np.where corresponds to the
    #   x-axis, while the second (<index> = 0) corresponds to the y-axis
    # NOTE: remember that pyfits/numpy invert x and y coordinates
    # NOTE: adding 1 to switch from array <index> to [pixel] coordinates

    x_min = ( np.min(np.where( image > 0 )[1]) + 1 )
    y_min = ( np.min(np.where( image > 0 )[0]) + 1 )
    x_MAX = ( np.max(np.where( image > 0 )[1]) + 1 )
    y_MAX = ( np.max(np.where( image > 0 )[0]) + 1 )
    

    # Checking if image can shall be restricted:
    
    if ADD_BUFFER is None:

      images_x_min.append(1)
      images_y_min.append(1)
      images_x_MAX.append(fits_NAXIS1[f])
      images_x_MAX.append(fits_NAXIS2[f])
    
    else:
    
      images_x_min.append(x_min) if (x_min > 1)              else images_x_min.append(1)
      images_y_min.append(y_min) if (y_min > 1)              else images_y_min.append(1)
      images_x_MAX.append(x_MAX) if (fits_NAXIS1[f] > x_MAX) else images_x_MAX.append(fits_NAXIS1[f])
      images_y_MAX.append(y_MAX) if (fits_NAXIS2[f] > y_MAX) else images_x_MAX.append(fits_NAXIS2[f])


  margin_r = max( [(images_x_MAX[f] - coords[f][0]   ) for f in range(len(fits_names))] )
  margin_l = max( [(coords[f][0]    - images_x_min[f]) for f in range(len(fits_names))] )
  margin_t = max( [(images_y_MAX[f] - coords[f][1]   ) for f in range(len(fits_names))] )
  margin_b = max( [(coords[f][1]    - images_y_min[f]) for f in range(len(fits_names))] )


  # Canvas size (+BUFFER):
  canvas_NAXIS1 = 2 * max(margin_r,margin_l) + 2 * BUFFER
  canvas_NAXIS2 = 2 * max(margin_t,margin_b) + 2 * BUFFER
  #-----------------------------------------------------------------------------

  # projecting each image on new canvases --------------------------------------
  # Shifts necessary to place desired pixel at the center of the new canvas:
  shift_dx = []
  shift_dy = []
  # [pixel]
  
  canvas_names = []

  for f, fits in enumerate(fits_names):

    shift_dx.append( int(round(canvas_NAXIS1/2. - coords[f][0])) )
    shift_dy.append( int(round(canvas_NAXIS2/2. - coords[f][1])) )
    # NOTE: By construction, canvas_NAXIS* are even, therefore shift_* will be
    #	    integers
    #	    Here the type cast is used to prevent a warning when sampling pixel
    #	    intervals from the images

    HDU_fits = pyfits.open(fits)
    image  = HDU_fits[fits_exts[f]].data
    header = HDU_fits[fits_exts[f]].header

    canvas_image = np.zeros((canvas_NAXIS2,canvas_NAXIS1))

    canvas_image[(images_y_min[f]+shift_dy[f]):(images_y_MAX[f]+shift_dy[f]),(images_x_min[f]+shift_dx[f]):(images_x_MAX[f]+shift_dx[f])] \
    = image[images_y_min[f]:images_y_MAX[f],images_x_min[f]:images_x_MAX[f]]
    # NOTE: now using array <index> (no need for +1) !

    canvas_names.append("canvas_"+fits)
    
    
    canvas_header = header
    canvas_header['CRPIX1'] = ( canvas_header['CRPIX1'] + shift_dx[f] )
    canvas_header['CRPIX2'] = ( canvas_header['CRPIX2'] + shift_dy[f] )


    # suppressing writeto warnings - - - - - - - - - - - - - - - - - - - - - - -
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pyfits.writeto(canvas_names[f], canvas_image, header=canvas_header, clobber=True)

  #-----------------------------------------------------------------------------
  
  return canvas_names
################################################################################


# Main #########################################################################
if __name__ == '__main__':

  # parsing input --------------------------------------------------------------
  program_name = sys.argv[0]
  # NOTICE: first entry of sys.argv is always the path to program name


  parser = argparse.ArgumentParser(description=
                                 'Expand the FOV of the input images in order to\
				 keep the desired pixel position of each input\
				 image at the center of the new canvas')

  parser.add_argument('-FILE', action='append', metavar='',dest='files', default=[],
                    help='add file to input list (add one entry for each file)')

  parser.add_argument('-CENTER', action='append', metavar='',dest='centers', default=[], nargs=2,
                    help='[alternative to -WCS]: "x y" of the desired center of each input image (add one entry for each file)')

  parser.add_argument('-WCS', action='append', metavar='',dest='WCS', default=[], nargs=2,
                    help='[alternative to -CENTER]: "RA Dec" of the desired center for all input images [deg]')

  parser.add_argument('-BUFFER', action='store', metavar='',dest='BUFFER', default=None,
                    help='[Optional]: desired minimal border around image (default=just use input images) [pixel]')

  parser.add_argument('--SILENT', action='store_true', default=False,
                    help='suppress STDOUT')

  args = parser.parse_args()
  #-----------------------------------------------------------------------------

  if not (args.SILENT):
    print "%s::" % program_name

  if not (args.WCS):
  # No inpout WCS specified
  
    centers = args.centers
    fit_on_canvas(args.files, centers, "IMAGE", BUFFER=args.BUFFER, SILENT=args.SILENT)
  
  else:
  # Input WCS specified

    centers = args.WCS
    fit_on_canvas(args.files, centers, "WCS", BUFFER=args.BUFFER, SILENT=args.SILENT)
  
################################################################################
