import matplotlib
import numpy as np

# the code currently (May 2020) gives a lot of matplotlib deprecation warnings
# which are probably cleaned up with matplotlibs version 3.3.
# until then we just turn off these warnings:
import warnings
warnings.warn("ignore")

# store original plot parameters so that we can revert:
ORIG_MATPLOTLIB_CONF = dict(matplotlib.rcParams)

def homogenise_plot(fig_width=None, fig_height=None, columns=1, fontsize=8):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this function before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}; 1 is default
    fontsze : 8 pt default
    """

    assert(columns in [1,2])

    if fig_width is None:
        # figure widths for one-column and two-column printing
        # (taken from https://www.elsevier.com/authors/author-schemas/artwork-and-media-instructions/artwork-sizing)
        if columns == 2:
            fig_width = 7
        else:
            fig_width = 3.5

    if fig_height is None:
        golden_mean = (np.sqrt(5.0) - 1.0) / 2.0    # Aesthetic ratio
        fig_height = fig_width * golden_mean # height in inches

    params = {'backend': 'ps',
              'text.latex.preamble':
              [ r'\usepackage{siunitx}',
                r'\usepackage[utf8]{inputenc}',
                r'\usepackage[T1]{fontenc}',
                r'\DeclareSIUnit \jansky {Jy}' ],
              'axes.labelsize' : fontsize,
              'axes.titlesize' : fontsize,
              'font.size': fontsize,
              'legend.fontsize' : fontsize,
              'xtick.labelsize' : fontsize,
              'ytick.labelsize' : fontsize,
              'axes.linewidth' : 1,
              'lines.linewidth' : 1,
              'text.usetex' : True,
              'figure.figsize' : [fig_width, fig_height],
              'font.family' : 'serif',
              'savefig.bbox' : 'tight',
              'savefig.dpi' : 300  # set to 600 for poster printing or PR
                                  # figures
    }

    matplotlib.rcParams.update(params)

def revert_params():
    """
    reverts any changes done to matplotlib parameters and restores
    the state before homogenize_plot was called
    """

    matplotlib.rcParams.update(ORIG_MATPLOTLIB_CONF)
