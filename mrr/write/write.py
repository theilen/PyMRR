

#Version history
#
#====================
#
#changes in 1.1
#-----------------
#- removed plot()
#- removed display()
#
#changes in 1.2
#-----------------
#- added adjust_window()
#   - added window_zero()
#   - added window_center()
#   - added set_mask()
#
#changes in 1.3
#----------------
#- added save_image()
#- added obsolete-warning to save_bitmap()
#
#changes in 1.3.1
#---------------
#- switched to relative imports

__version__ = '1.3.1'
# $Source$


import numpy as np
import scipy.misc
import matplotlib.pyplot as plt


from .. import mrrcore as core
from ..plotting import display_plain


def adjust_window(images, window='center', mask='center'):
    '''
    Adjust the greyvalues of all given images, so that they all cover the same
    range of values.
    
    window : ['center' | 'zero']
        the method used to adjust the values:
            'center' : the center-value of all adjusted images is the same.
            'zero' : all images range from zero to the maximal value of all 
                     adjusted images.
                     
    mask : ['center' | 'min' | 'max' | float]
        the value assigned to masked pixels:
            'center' : exactly between the minimum and the maximum value of 
                       all images after adjustment.
            'min' : the minimal value of all images after adjustment.
            'max' : the maximal value of all images after adjustment.
            float : a specific value. Note that decimal-numbers use '.' 
                    instead of ','!
    '''
    if not window in ['zero', 'center']:
        raise AttributeError('No window "%s" defined!' % window)
    single = False
    gmin = 999.
    gmax = -999.
    grange = 0.
    if hasattr(images, 'shape'):
            if len(images.shape) == 2: #single image
                single = True
                gmin = core.min(images)
                gmax = core.max(images)
                grange = gmax - gmin
            elif len(images.shape) != 3:
                raise TypeError('Invalid dimensions for image data')
    if not single:
        for item in images:
            gmin = min(gmin, core.min(item))
            gmax = max(gmax, core.max(item))
            grange = max(grange, core.max(item)-core.min(item))
    print 'adjust_window: gmin,gmax,grange =',gmin, gmax, grange
    if single:
        return eval('window_' + window
                    + '(images,%f,%f,%f' % (gmin,gmax,grange)        
                    + ',"%s")' % str(mask)
                    )
    else:
        result = []
        for item in images:
            result.append( eval('window_' + window
                                + '(item,%f,%f,%f' % (gmin,gmax,grange)
                                + ',"%s")' % str(mask)
                                ) 
                          )
        return result

    
def set_mask(array, mask, valrange):
    if mask=='center':
        value = valrange/2.
    elif mask=='min':
        value = 0.
    elif mask=='max':
        value = valrange
    else:
        value = float(mask)
    array['phase'][array['mask']==False] = value
    
def window_zero(array, minval, maxval, valrange, mask):
    item = core.sub(array, core.min(array))
    set_mask(item, mask, valrange)
    item['phase'][:,0] = valrange
    return item
    
def window_center(array, minval, maxval, valrange, mask):
    min_ = core.min(array)
    max_ = core.max(array)
    space = (valrange - (max_ - min_))/2.
    corr = space - min_
    item = core.add(array, corr)
    set_mask(item, mask, valrange)
    item['phase'][:,0] = valrange
    item['phase'][:,1] = 0.
    return item
        
    
    

def save_txt(path, array, field='phase'):
    '''
    Saves the values in a specified field of array as a matrix in a textfile.
    Only meaningful with 1- or 2-dimensional data!

    The created file with override existing files without notice.

    Parameters
    ----------
    path : str
        String containing the path to the desired filename. An ending is added
        automatically.
    array : MRRArray
        array containing the data. Only works with one or two dimensional data.
    field : str
        default: 'phase'. Field to save to file.
    '''
    try:
        np.savetxt(path, array[field])
    except ValueError:
        raise ValueError('Array does not contain a field %s' % field)


def create_plotdata(*plotdata):
    pdata = []
    edata = []
    labels = []
    for data_set in plotdata:
        try:
            if type(data_set[1]) == type(''):
                label = data_set[1]
                data = data_set[0]
            else:
                label = ''
                data = data_set
        except IndexError:
            label = ''
            data = data_set
        pdata.append(data['phase'])
        edata.append(data['dev'])
        labels.append(label)
    return labels, pdata, edata



 
def plot_to_txt(filename, *plot_data, **keyargs):
    '''
    save plot-data to a textfile in gnuplot style.
    
    plot_data : sequence
        lists or 1-dimensional arrays containing the data to be saved. Every
        list contains the data of one column. Make sure all provided data-lists
        conatin the same amount of elements! 
    filename : str
        relative or absolute path to the file the data is saved to.
    docstring : str (optional)
        Written to the first lines of the resulting file. Every line will be 
        marked as commentary by a preceeding '#'
    
    **Tested**
    '''
    labels, pdata, edata, = create_plotdata(*plot_data)
    with open(filename, 'w') as outfile:
        #optinal docstring
        try:
            docstring = keyargs['docstring']
        except KeyError:
            pass
        else:            
            docstring = '#' + docstring.replace('\n','\n#') + '\n'
            outfile.write(docstring)
        #label
        line = '#'
        for label in labels:
            if label == '':
                label = 'n/a'
            line += '%s\t' % label
        outfile.write(line + '\n')
        #data
        for i in range(len(pdata[0])):
            line = ''
            for d in range(len(pdata)):
                line += '%e\t%e\t' % (pdata[d][i], edata[d][i])
            outfile.write(line + '\n')
        outfile.write('\n#EOF')                
        
        
def save_bitmap(filename, array, field='phase'):
    '''
    Saves the array a 2-d-array as bitmap.
    
    The resulting image will not contain the original precision of the data!
    
    Corel wird beim importieren der erzeugten Datei vermutlich die Groesse 
    nicht richtig verarbeiten. Evtl hilft es dann, die Datei mit paint zu
    oeffnen und neu zu speichern.
    
    ACHTUNG: Aus den vorgenannten Gruenden sollte stattdessen die Funktion 
        save_image (in mrr.write)
    genutzt werden!
    '''
    scipy.misc.imsave(filename, array[field])
    print 'This function is obsolete and only included for backward-compatibility! Use save_image instead!'
    return filename
    

def save_image(filename, img, field='phase', figwidth=4.0,
               crop=False, crop_range=10,
               **kwargs):
    '''
    Saves a 2D-Array to an image-file.
    
    The file-format is derived automatically from filename.
    Allowed formats are, a.o.: .pdf and .png.

    ACHTUNG: Bitmap-Dateien (.bmp) koennen mit dieser Funktion NICHT erstellt 
    werden! png ist aber in jeder Hinsicht gleichwertig.
    
    
    img : MRRArray
        img-data to be visualized
    field : ['phase' | 'dev' | 'mask'] (optional)
        Which field to plot. Default: 'phase'
    figwidth : skalar (optional)
        Width of the resulting image in inch. Default: 4.0
        
    crop : [True | False] (optional)
        Wether to crop the image based on a['mask']
    crop_range : int (optional)
        Adjusts the space around mask to be shown as well. Default: 10
        
    The plotting of the image is done by pyplot.imshow(). This can be 
    custumized through additional key-arguments, that are passed through to 
    imshow(). In particular, this can be used to adjust the colormap:
    
    cmap=mycolormap uses your own custom colormap. available colormaps can be 
    found online.
    Using vmin and vmax you can adjust the minimum and maximum values of the 
    colormap, therefore adjusting the windowing of the image.
    '''
    plt.interactive(False)

    #fig = plt.figure()
    display_plain(img, field=field, figwidth=figwidth,
                  crop=crop, crop_range=crop_range,
                  **kwargs)
    plt.savefig(filename)
    
    plt.close(plt.gcf())
    plt.interactive(True)
