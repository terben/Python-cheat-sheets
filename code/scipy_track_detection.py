#!/usr/bin/env python3

# This script identifies satellite tracks in astronomical
# images as discussed in the corresponding notebook.
# If one or more tracks are found, a plot is created from
# the original image and the detected tracks.

# necessary modules for this sheet:
import numpy as np

# to perform object detection and analysis
import scipy.ndimage.measurements as snm
import scipy.ndimage.filters as snf
import scipy.ndimage.morphology as snmo
import scipy.linalg as sl

# to handle FITS image I/O and robust image statistics
import astropy.io.fits as aif
import astropy.stats as ast

# to manage plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mp

# settings for nicer plots:
# font size of labels etc,
matplotlib.rcParams['font.size'] = 18
# line width of coordinate axes
matplotlib.rcParams['axes.linewidth'] = 2.0

# FREE PARAMETERS YOU CAN SET:
image = '../data/track_demo.fits'  # the image to analyse - put an own image
                                # here if you like
out_plot = 'track_demo.png'     # name of teh output plot
k = 1.0                         # we detect objects having pixels
                                # with k * sigma_sky
req_pix = 30                    # The number of pixels that must be above
                                # k * sigma_sky to be included in the
                                # analysis

# load image
hdu = aif.open(image)
data = hdu[0].data

# obtain 'robust statistics' around the sky-background of the image.
# They is needed for image plotting below
data_mean, data_median, data_sigma = ast.sigma_clipped_stats(data[::10,::10])
min_val = data_median - 3.0 * data_sigma
max_val = data_median + 10.0 * data_sigma

print("statistics from %s" % (image))
print("mean: %.2f; median: %.2f; sigma: %.2f" %
      (data_mean, data_median, data_sigma))

# smooth image for object detection:

# pixels exceeding 1 \sigma_sky without any smoothing:
sigma_gauss = 1.0 # smoothing sigma
data_smooth = snf.gaussian_filter(data, sigma=sigma_gauss)
det_smooth = data_smooth > k * data_sigma

# label connected features in the detection image
objects, n_objects = snm.label(det_smooth)

print("%d candidate objects found" % (n_objects))

# obtain slices from the connected features (sources) identified with the labels
# and analyse everything having 30 pixels or more
slices = snm.find_objects(objects)

track_pixels = np.zeros(objects.shape, bool)

# now check whether there are tracks among the objects:
for i, slc in enumerate(slices):
    # Note that we only work with the detection image, not on the original
    # pixel data here!
    det_pixels = objects[slc]

    # cehck whether the object has a significant size -
    # a necessary condition for a track:
    slice_height, slice_width = det_pixels.shape

    # The '300' is the width/height of the image
    if slice_height > 300 / np.sqrt(2.) or slice_width > 300 / np.sqrt(2.):
        det_pixels = (det_pixels == i + 1)
        x = np.arange(slice_width)
        y = np.arange(slice_height)
        flux = np.sum(det_pixels)
        x_cen = np.sum(x[np.newaxis,:] * det_pixels) / flux
        y_cen = np.sum(y[:,np.newaxis] * det_pixels) / flux
        c_xx = np.sum(((x - x_cen)**2)[np.newaxis,:] * det_pixels) / flux
        c_yy = np.sum(((y - y_cen)**2)[:,np.newaxis] * det_pixels) / flux
        c_xy = np.sum((x - x_cen)[np.newaxis,:] *
                      (y - y_cen)[:,np.newaxis] * det_pixels) / flux
        C = np.array([[c_xx, c_xy], [c_xy, c_yy]])
        w, v = sl.eig(C)
        a, b = np.max(2.0 * np.sqrt(np.real(w))), \
               np.min(2.0 * np.sqrt(np.real(w)))

        # The final check whether we have a track:
        if a > (300 / (2. * np.sqrt(2.))) and b / a < 0.1:
            track_pixels[slc] = det_pixels

# make a plot from original image and track pixels if track(s) were found:
if np.sum(track_pixels) > 0:
    print('Track(s) detected! Creating plot %s' % (out_plot))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    ax1.imshow(data, origin='lower', cmap='gray_r',
              vmin=min_val, vmax=max_val)
    ax1.set_xlabel('x [pix]')
    ax1.set_ylabel('y [pix]')
    ax1.set_title('original data')
    ax2.imshow(data, origin='lower', cmap='gray_r',
              vmin=min_val, vmax=max_val)
    # dilate the track by 5 on its borders pixels to be sure to catch all
    # pixels belonging to it
    track_image = np.zeros(track_pixels.shape + (4,), np.uint8)
    newtrack_pixels = snmo.binary_dilation(track_pixels,
                                           structure=np.ones((5,5)))
    track_image[:,:,2] = (newtrack_pixels > 0) * 255
    track_image[:,:,3] = (newtrack_pixels > 0) * 255
    ax2.imshow(track_image, origin='lower')
    ax2.set_xlabel('x [pix]')
    ax2.set_ylabel('y [pix]')
    ax2.set_title('detected track')

    plt.tight_layout()
    plt.savefig(out_plot)
else:
    print('No track in %s dound.' % (image))
