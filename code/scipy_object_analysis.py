#!/usr/bin/env python3

# This script analyses objects in an optical astronomical FITS image
# as discussed in the accompanying notebook.
# see below for some free parameters that you can adjust for own experiments.

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
k = 3.0                         # we detect objects having pixels
                                # with k * sigma_sky
req_pix = 10                    # The number of pixels that must be above
                                # k * sigma_sky to be included in the
                                # analysis

#image = '../tracks/track_demo.fits'
#k = 5.0
#req_pix = 10.0

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

# plot image and smoothed detection image:
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
ax1.imshow(data, origin='lower', cmap='gray_r',
           vmin=min_val, vmax=max_val)
ax1.set_xlabel('x [pix]')
ax1.set_ylabel('y [pix]')
ax2.imshow(det_smooth, origin='lower', cmap='gray')
ax2.set_xlabel('x [pix]')
ax2.set_ylabel('y [pix]')
plt.tight_layout()
plt.savefig('data_det.png', dpi=200)
plt.close(fig)

# label connected features in the detection image
objects, n_objects = snm.label(det_smooth)

print("%d candidate objects found" % (n_objects))

# obtain slices from the connected features (sources) identified with the labels
# and analyse everything having 30 pixels or more
slices = snm.find_objects(objects)

for i, curr_slice in enumerate(slices):
    curr_object = objects[curr_slice]
    curr_object = (curr_object == i + 1)

    n_pix = np.sum(curr_object)

    if n_pix > req_pix:
        print("Analysing object %d with %d pixels of %d sigma above the sky." %
              (i, n_pix, k))
        obj_det = curr_object
        obj_data = data[curr_slice]
        obj_pixels = obj_data * obj_det

        slice_height, slice_width = obj_det.shape
        x = np.arange(slice_width)
        y = np.arange(slice_height)

        # obtain object quantities:
        flux = np.sum(obj_pixels)

        # estimate object center:
        x_cen = np.sum(x[np.newaxis,:] * obj_pixels) / flux
        y_cen = np.sum(y[:,np.newaxis] * obj_pixels) / flux

        # estimate second order moments:
        c_xx = np.sum(((x - x_cen)**2)[np.newaxis,:] * obj_pixels) / flux
        c_yy = np.sum(((y - y_cen)**2)[:,np.newaxis] * obj_pixels) / flux
        c_xy = np.sum((x - x_cen)[np.newaxis,:] *
                      (y - y_cen)[:,np.newaxis] * obj_pixels) / flux

        # We assume that our object can be well represented by an ellipse
        # and we estimate the ellipse axes and its position angle:
        C = np.array([[c_xx, c_xy], [c_xy, c_yy]])

        # obtain the eigenvalues and eigenvectors of C and extract
        # desires ellipse quantities:
        w, v = sl.eig(C) # w: complex eigenvalues; v: eigenvectors

        # ellipse axes:
        a, b = 2.0 * np.sqrt(np.real(w))
        a, b = np.max(2.0 * np.sqrt(np.real(w))), \
               np.min(2.0 * np.sqrt(np.real(w)))

        # ellipse position angle:
        phi = np.rad2deg(np.arctan2(v[1][0], v[0][0]))
        phi = np.rad2deg(0.5 * np.arctan2(2.0 * c_xy, c_xx - c_yy))

        # make a plot of the object and its ellipse:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        ax.imshow(data[curr_slice], origin='lower', cmap='gray_r',
                  vmin=min_val, vmax=max_val)
        ellipse = mp.Ellipse((x_cen, y_cen),
                             2.0 * a, 2.0 * b, angle=phi, linewidth=5,
                             edgecolor='lightgreen', fill=False)
        ax.add_patch(ellipse)
        ax.set_xlabel('x [pix]')
        ax.set_ylabel('y [pix]')
        plt.tight_layout()
        plt.savefig("object_%d.png" % (i), dpi=200)
        plt.close(fig)
