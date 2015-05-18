# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:22:36 2014

last change: Sun Nov 09 2014

@author: Sebastian Theilenberg
"""

__version__ = '1.4'
# $Source$

import numpy as np
import matplotlib.pyplot as plt


_default_clabel_tex = {
    "mask": r"mask",
    "phase": r"$\phi\,\,\left[\frac{\mathrm{rad}}{2\pi}\right]$",
    "dev": r"$\Delta\phi\,\,\left[\frac{\mathrm{rad}}{2\pi}\right]$"
    }

_default_clabel_wtex = {
    "mask": "mask",
    "phase": "phase [rad/2pi]",
    "dev": "std. dev. phase [rad/2pi]"
    }


def _cbar_label(field):
    try:
        if plt.rcParams["text.usetex"] is True:
            string = _default_clabel_tex[field]
        else:
            string = _default_clabel_wtex[field]
    except KeyError:
        string = ""
    finally:
        return string


def display(img, field='phase',
            crop=False, crop_range=10,
            shrink_cbar=.9,
            grid=True,
            hold=False,
            ckwargs=None,
            clabel=None,
            **kwargs):
    '''
    Display data of an MRRArray.

    If dim(img) = 3 the first data-set in axis 0 is used.
    On default the phase-data is displayed, other fields may be displayed by
    providing <field>. The image may be zoomed in based on img['mask']. To do
    so, set crop=True.

    The actual plotting of the data is done using plt.imshow(). All additional
    key-arguments are passed through, so if you want to modify e.g. the
    colormap, please see the documentation of pyplot.

    Parameters
    ----------
    img : MRRArray
        img-data to be visualized
    field : ['phase' | 'dev' | 'mask'] (optional)
        Which field to plot. Default: 'phase'

    hold : [True | False] (optional)
        If True, the figure used to image the data will stay active, so that
        additional data can be plotted to the same image.
        Default: False

    crop : [True | False] (optional)
        Wether to crop the image based on a['mask']
    crop_range : int (optional)
        Adjusts the space around mask to be shown as well. Default: 10
    shrink_cbar : value (optional)
        Adjusts the size of the colorbar
    grid : [True | False] (optional)
        Whether to draw a grid
    ckwargs : dictionary
        optional arguments to be passed to the colorbar() command.
    clabel : string (optional)
        Label for the colorbar. Set to "" to ignore labeling.

    Returns
    -------
    fig : matplotlib.figure
        Only if hold==True. The figure-object containing the image.
    cbar : matplotlib.colorbar
    '''
    fig = plt.figure()
    # crop 3d-arrays
    if len(img.shape) == 3:
        img = img[0]
    # crop image based on mask
    if crop:
        try:
            ylims, xlims = crop_image(img, crop_range)
        except ValueError:
            print 'Cropping not possible. Is there a mask present?'
        else:
            plt.xlim(xlims)
            plt.ylim(ylims[::-1])

    # get phase-data
    try:
        temp = img[field].view(np.ndarray)
    except ValueError:
        print 'Field \'%s\' not found!' % field
        temp = img

    # plot
    imshow_dict = {'cmap': plt.get_cmap('Greys_r'),  # colormap
                   'interpolation': 'none'           # interpolation
                   }
    imshow_dict.update(kwargs)
    plt.imshow(temp, **imshow_dict)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='out')

    # colorbar
    colorbar_dict = {'shrink': shrink_cbar}
    if ckwargs is not None:
        colorbar_dict.update(ckwargs)
    cbar = plt.colorbar(shrink=shrink_cbar)
    if clabel is None:
        clabel = _cbar_label(field)
    if clabel:
        cbar.set_label(clabel)

    if grid:
        plt.grid()

    if not hold:
        plt.show()  # setzt figure zurück
    else:
        return fig, cbar


def display_plain(img, field='phase',
                  crop=False, crop_range=10,
                  figwidth=4.0,
                  grid=True,
                  hold=True,
                  **kwargs):
    '''
    Display data of an MRRArray without additional elements, e.g. without axes,
    colorbar, etc. This creates an output as the one of save_image() without
    saving the image.

    All meaningful arguments work like the ones of display().

    img : MRRArray
        img-data to be visualized
    field : ['phase' | 'dev' | 'mask'] (optional)
        Which field to plot. Default: 'phase'

    hold : [True | False] (optional)
        If True, the figure used to image the data will stay active, so that
        additional data can be plotted to the same image.
        Default: True

    crop : [True | False] (optional)
        Wether to crop the image based on a['mask']
    crop_range : int (optional)
        Adjusts the space around mask to be shown as well. Default: 10

    The actual plotting of the data is done using plt.imshow(). All additional
    key-arguments are passed through, so if you want to modify e.g. the
    colormap, please see the documentation of pyplot.
    '''
    imshow_dict = dict(interpolation='none',
                       cmap='Greys_r')
    imshow_dict.update(**kwargs)

    # crop 3d-arrays
    if len(img.shape) == 3:
        img = img[0]

    try:
        temp = img[field].view(np.ndarray)
    except ValueError:
        print 'Field \'%s\' not found!' % field
        temp = img

    fig = plt.figure()

    if crop:
        try:
            ylims, xlims = crop_image(img, crop_range)
        except ValueError:
            print 'Cropping not possible. Is there a mask present?'
            crop = False
    if crop is False:
        y_max, x_max = img.shape
        xlims, ylims = (0, x_max-1), (0, y_max-1)
    plt.xlim(xlims)
    plt.ylim(ylims[::-1])

    ratio = calc_ratio(ylims, xlims)

    fig.set_size_inches(figwidth, figwidth*ratio)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims[::-1])
    fig.add_axes(ax)

    ax.imshow(temp, **imshow_dict)

    if not hold:
        plt.show(fig)


def crop_image(array, crop_range=10):
    '''
    Calculates x- and y-limits based on the masked of the provided array.

    Returns ylims, xlims (tuple).
    If there is no masked present, ValueError is raised.
    '''
    if not valid_mask(array):
        raise ValueError("No masked present")

    indices = unmasked_indices(array)
    ylims, xlims = crop_limits(indices, crop_range)

    return ylims, xlims


def overview(array, field="phase", nx=5, imageaxis=0, averageaxis=0,
             autoadjust=False):
    """
    Creates a large overview plot displaying all images in array in a grid.
    Written to be used with the data returned from read_dicom_set.

    parameters:
    -----------
    array : array-like
        data to be displayed. This may be a multi-dimensional array (at least
        3D) or a list of arrays.
    field : str, optional, default: "phase"
        field to be displayed when using MRRArrays
    nx : int, optional, default: 5
        number of columns
    imageaxis : int, optional, default: 0
        axis of array to iterate over to find images
    averageaxis : int or None, optional, default: 0
        axis along to average the elements of array. If averaging is not
        wanted, set to None
    autoadjust : bool, optional, default: False
        If set to True, all images span the same range of grey values.
    """
    data = np.asarray(array)
    # roll image axis to front
    np.rollaxis(data, imageaxis, 0)

    # determine whether MRRArrays are present
    try:
        data[0][field]
    except IndexError:
        field = None

    if autoadjust:
        if field is not None:
            vmin = min([d_.dataview(field).min() for d_ in data])
            vmax = max([d_.dataview(field).max() for d_ in data])
        else:
            vmin = min([d_.min() for d_ in data])
            vmax = max([d_.max() for d_ in data])
    else:
        vmin = vmax = None

    n = data.shape[0]
    ny = int(np.ceil(1.0*n/nx))
    fig, axarr = plt.subplots(ny, nx, figsize=(nx*4.0, ny*3.5))

    i = 0
    for y in xrange(ny):
        for x in xrange(nx):
            ax = axarr[y, x]
            # create average
            try:
                img = data[i]
            except IndexError:
                ax.axis("off")
                continue
            if averageaxis is not None:
                img = img.mean(axis=averageaxis)
            else:
                if img.dims >= 3:
                    img = img[0]
            # extract data as numpy array
            if field is not None:
                img = img[field].view(np.ndarray)

            ax.imshow(img, interpolation='None', cmap='Greys_r',
                      vmin=vmin, vmax=vmax)
            ax.set_title("Image {}".format(i))

            i += 1


def parse_pd(pd):
    '''
    Parse an item of ( [x|(x0,x1),] y [, 'fmt'] ) and return the three possible
    fields:

    x : [None | (x_min, x_max) | seq of x-values]
    y : sequence of data
    fmt : [None | format string]
    '''
    item = pd.pop()
    try:
        fmt = '' + item
    except TypeError:
        fmt = ''
    else:
        item = pd.pop()
    finally:
        y = item
        try:
            x = pd.pop()
        except IndexError:
            x = None
    return x, y, fmt


def valid_mask(array):
    '''
    Returns True if array has a field 'mask' and at least one element has this
    field set to 'False', and False otherwise.
    '''
    try:
        if np.any(array.mask == False):
            return True
    except ValueError:
        pass
    return False


def unmasked_indices(array):
    '''
    Returns a tuple of indizes of y, where y['mask'] == True:

        ([[...,] yindices,] xindices,)
    '''
    indices = np.where(array.mask == True)

    return indices


def crop_limits(indices, border):
    'Returns (y_min,y_max),(x_min,x_max) values based on unmasked_indices'
    y_min = max(0, indices[0].min() - border)
    y_max = indices[0].max() + border
    x_min = max(0, indices[1].min() - border)
    x_max = indices[1].max() + border
    return (y_min, y_max), (x_min, x_max)


def calc_ratio(ylims, xlims):
    'calculates y/x'
    y = ylims[1] - ylims[0]
    x = xlims[1] - xlims[0]
    return 1.0*y/x


def plot(*args, **keyargs):
    """
    Plot lines of MRRArrays in a pyplot-fashion.

    *args* is of variable length containing *x*, *y* pairs with an optional 
    format string 'fmt'. Each of the following is legal:

        plot(y)                 #plot y with x being the original index
        plot(y, 'g-')           #ditto with a green line

        plot((0,50), y)         #plot y from 0 to 50 (x being the orig. index)
        plot ((0,50), y, 'bo')  #ditto with blue circle markers
        
        plot(x, y)              #plot y versus x
        plot(x, y, 'r+')        #ditto with red plusses
    
    An arbitrary number of y,'fmt' groups can be specified as in:
    
        plot(y1, 'go', y2, 'r+')
        
    Every y will be plotted using y['phase'] as value and y['dev'] for 
    errorbars. By default, datapoints with y['mask']==False are ignored.
    
    The actual plotting is done with matplotlib.pyplot.errorbar(). The 
    formatstring and all keyword arguments work here just like there. 
    
    Note that keywordarguments apply to all lines made with this plot command,
    e.g.:
    
        plot(y1, 'r+', y2, '', label='lines')
        
    will assign the label 'lines' to both y1 and y2.
    
    Additionally the following keyword arguments are supported:
    
        mask_pixels : bool
            If set to False, all pixels will be plotted, regardless of 
            y['mask']. Is a sequence of x-values is provided, this value is 
            set to False. (Default: False)
        
        errors: bool
            If set to False, no errorbars are plotted. (Default: True)
            
        labels : bool
            Whether to label the axes. (Default: False)
            
        return_index : bool
            if set to True, the values x_min, x_max are returned specifying 
            the plotting limits determined. Useful if further plotting needs 
            to be align with the phase data. (Default: False)
    """
    #handle positional arguments - sort plot-data
    if len(args) > 3:
        plots = [list(args[i:i+2]) for i in range(0, len(args), 2)]
    else:
        plots = [list(args)]
    #handle custom keyarguments
    try:
        mask_pixels = keyargs.pop('mask_pixels')
    except KeyError:
        mask_pixels = True
    try:
        errors = keyargs.pop('errors')
    except KeyError:
        errors = True
    try:
        return_index = keyargs.pop('return_index')
    except KeyError:
        return_index = False
    try:
        labels = keyargs.pop('labels')
    except KeyError:
        labels = False

    #handle plot-data
    for pd in plots:
        x, y, fmt = parse_pd(pd)
        #create indices
        if not x==None and not isinstance(x, tuple): #array of indices provided
            #create index-array            
            indices = np.arange(0, len(x), dtype=np.int)
            mask_pixels = False
            if indices.shape != y.shape:
                print 'Data does not match amount of x-values! Check your values!'
            #x-values for plotting
            xvals = x
        else:
            try: #Min,Max provided
                x_min, x_max = x
            except TypeError: #no x provided
                x_min = 0
                x_max = len(y)
            finally:#create index-array
                indices = np.arange(x_min, x_max, dtype=np.int)
            #delete masked data
            if mask_pixels:
                masked_indices, = unmasked_indices(y)
                indices = np.asarray(
                                [x for x in indices if x in masked_indices],
                                dtype=np.int)
            #x-values for plotting (=indices)
            xvals = indices.copy()
        #plot
        values = y.phase[indices]
        if errors:
            err = y.dev[indices]
        else:
            err = None
        plt.errorbar(xvals, values, yerr=err, fmt=fmt, **keyargs)

        #labels
        if labels:
            if plt.rcParams['text.usetex'] == True:
                label = _plot_labels['tex']
            else:
                label = _plot_labels['norm']
            plt.ylabel(label['yl'])
            plt.xlabel(label['xl'])

        if return_index:
            return indices[0], indices[-1]


_plot_labels = {'tex' : {'yl' : r'$\phi\,\left[rad/2\pi}\right]$',
                         'xl' : r'Pixel'},
                'norm' : {'yl' : r'phase [rad/2pi]',
                          'xl' : r'Pixel'}
                }
