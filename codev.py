"""
Author:         David Henry Lippman
File:           codev.py
Date created:   07/21/20
Date modified:  07/21/20

"""

import numpy as np
from win32com.client import DispatchWithEvents
import pythoncom
import sys


class ICVApplicationEvents:

    def OnLicenseError(self, error):
        # This event handler is called when a licensing error is
        # detected in the CODE V application.
        print("License error: %s " % error)

    def OnCodeVError(self, error):
        # This event handler is called when a CODE V error message is issued
        print("CODE V error: %s " % error)

    def OnCodeVWarning(self, warning):
        # This event handler is called when a CODE V warning message is issued
        print("CODE V warning: %s " % warning)

    def OnPlotReady(self, filename, plotwindow):
        # This event handler is called when a plot file, refered to by filename,
        # is ready to be displayed.
        # The event handler is responsible for saving/copying the
        # plot data out of the file specified by filename
        print("CODE V Plot: %s in plot window %d" % (filename, plotwindow))


def init_cv(cd="C:\\CVUSER"):

    """


    Initializes a CODE V COM session.


    """

    # Try connecting

    try:
        cv = DispatchWithEvents("CodeV.Application", ICVApplicationEvents)
        cv.StartingDirectory = cd
        cv.StartCodeV()
        return cv

    # Otherwise show error

    except pythoncom.com_error as error:
        print("---------")
        args = error.args
        for a in args:
            print(a)
        print(error.strerror)
        stop_cv(cv)


def stop_cv(cv):

    """


    Closes a CODE V COM session.


    """

    cv.StopCodeV()
    del cv


def cmd(s, out=True):
    
    """
    
    
    Issues CODE V command.
    
    
    s:          command string

    out:        flag for printing output
    
    """

    print(s)
    result = cv.Command(s)
    if out:
        print(result)
    return result


def eva(s, dtype=float):
    
    """
    
    
    Evaluates CODE V database item.
    
    
    s:          evaluation string without database open-close parentheses

    dtype:      data type of output to return
    
    """

    return dtype(cv.EvaluateExpression('(' + s + ')'))


def vec_to_str(vec, delim=' '):

    """


    Returns a string of the numerical contents of a vector, separated by some
    delimiter.


    vec:        (m,) vector to return the numerical contents of

    delim:      delimiter; default is single space

    """

    vec_str = ''

    for v in vec:
        vec_str += delim + '{0:0.12f}'.format(v)

    return vec_str


def set_num_z(num_z):

    """


    Ensures there are (at least) the specified number of zoom positions. If not,
    adds additional zoom positions.


    num_z:          desired number of zoom positions

    """

    while num_z > eva('NUM Z'):
        cmd('INS Z{0:0.0f}+1'.format(eva('NUM Z')))

    return


def set_num_f(num_f):

    """


    Ensures there are (at least) the specified number of field positions. If
    not, adds additional field positions.


    num_zoom:       desired number of zoom positions

    """

    while num_f > eva('NUM F'):
        cmd('INS F{0:0.0f}+1'.format(eva('NUM F')))

    return


def set_pupil(epd=None, fno=None, na=None, nao=None):
    
    """
    
    
    Sets the system pupil settings.
    
    
    epd:        entrance pupil diameter [mm]; scalar (all zoom positions) or
                array (zooms parameter)

    fno:        f/number; scalar (all zoom positions) or array (zooms parameter)

    na:         image space numerical aperture; scalar (all zoom positions) or
                array (zooms parameter)

    nao:        object space numerical aperture; scalar (all zoom positions) or
                array (zooms parameter)
     
    """

    # Set EPD, if provided

    if epd is not None:

        # Update EPD for all zoom positions

        if np.isscalar(epd):
            cmd('EPD {0:0.12f}'.format(epd))

        # Update EPD for different zoom positions

        else:
            set_num_z(epd.size)
            cmd('ZOO EPD')
            for z, val in enumerate(epd):
                cmd('EPD Z{0:0.0f} {1:0.12f}'.format(z + 1, val))

    # Set f/number, if provided

    if fno is not None:

        # Update f/number for all zoom positions

        if np.isscalar(fno):
            cmd('FNO {0:0.12f}'.format(fno))

        # Update f/number for different zoom positions

        else:
            set_num_z(fno.size)
            cmd('ZOO FNO')
            for z, val in enumerate(fno):
                cmd('FNO Z{0:0.0f} {1:0.12f}'.format(z + 1, val))

    # Set NA, if provided

    if na is not None:

        # Update NA for all zoom positions

        if np.isscalar(na):
            cmd('NA {0:0.12f}'.format(na))

        # Update NA for different zoom positions

        else:
            set_num_z(na.size)
            cmd('ZOO NA')
            for z, val in enumerate(na):
                cmd('NA Z{0:0.0f} {1:0.12f}'.format(z + 1, val))

    # Set object space NA, if provided

    if nao is not None:

        # Update object space NA for all zoom positions

        if np.isscalar(nao):
            cmd('NAO {0:0.12f}'.format(nao))

        # Update object space NA for different zoom positions

        else:
            set_num_z(nao.size)
            cmd('ZOO NAO')
            for z, val in enumerate(nao):
                cmd('NAO Z{0:0.0f} {1:0.12f}'.format(z + 1, val))


def set_wl(wl=None, ref=1):
    
    """
    
    
    Sets the system wavelength settings


    wl:         (m,) array of wavelengths to set [nm]

    ref:        reference wavelength integer (zero-indexed for python)

    """

    # Update wavelengths, if provided

    if wl is not None:

        # Ensure wavelengths are in a numpy array

        wl = np.array(wl)

        # Sort (descending) the wavelength(s)

        wl = np.sort(wl)[::-1]

        # Issue wavelength command

        cmd('WL' + vec_to_str(wl))

    # Set reference wavelength

    cmd('REF {0:0.0f}'.format(ref + 1))


def set_fld(yan=None, xan=None, yim=None, xim=None, yri=None, xri=None,
            yob=None, xob=None, set_vig=True):

    """


    Sets the system field-of-view settings.


    yan:        Y object angle [deg]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    xan:        X object angle [deg]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    yim:        Y paraxial image height [l.u.]; vector (all zoom positions) or
                array (zooms parameter) where rows correspond to different zooms

    xim:        X paraxial image height [l.u.]; vector (all zoom positions) or
                array (zooms parameter) where rows correspond to different zooms

    yri:        Y real image height [l.u.]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    xri:        X real image height [l.u.]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    yob:        Y object height [l.u.]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    xob:        X object height [l.u.]; vector (all zoom positions) or array
                (zooms parameter) where rows correspond to different zooms

    set_vig:    flag for setting vignetting values. Uses CV_MACRO:setvig.seq
                macro.

    """

    # Sets Y object angle, if provided

    if yan is not None:

        # Update Y object angle for all zoom positions

        if yan.ndim == 1:
            cmd('YAN' + vec_to_str(yan))

        # Update Y object angle for different zoom positions

        elif yan.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(yan.shape[1])
            set_num_z(yan.shape[0])

            # Loop over field

            for f in range(yan.shape[1]):

                # Zoom field position

                cmd('ZOO YAN F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(yan.shape[0]):

                    # Set field value

                    cmd('YAN F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   yan[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets X object angle, if provided

    if xan is not None:

        # Update X object angle for all zoom positions

        if xan.ndim == 1:
            cmd('XAN' + vec_to_str(xan))

        # Update X object angle for different zoom positions

        elif xan.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(xan.shape[1])
            set_num_z(xan.shape[0])

            # Loop over field

            for f in range(xan.shape[1]):

                # Zoom field position

                cmd('ZOO XAN F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(xan.shape[0]):
                    # Set field value

                    cmd('XAN F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   xan[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets Y paraxial image height, if provided

    if yim is not None:

        # Update Y paraxial image height for all zoom positions

        if yim.ndim == 1:
            cmd('YIM' + vec_to_str(yim))

        # Update Y paraxial image height for different zoom positions

        elif yim.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(yim.shape[1])
            set_num_z(yim.shape[0])

            # Loop over field

            for f in range(yim.shape[1]):

                # Zoom field position

                cmd('ZOO YIM F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(yim.shape[0]):
                    # Set field value

                    cmd('YIM F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   yim[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets X paraxial image height, if provided

    if xim is not None:

        # Update X paraxial image height for all zoom positions

        if xim.ndim == 1:
            cmd('XIM' + vec_to_str(xim))

        # Update X paraxial image height for different zoom positions

        elif xim.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(xim.shape[1])
            set_num_z(xim.shape[0])

            # Loop over field

            for f in range(xim.shape[1]):

                # Zoom field position

                cmd('ZOO XIM F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(xim.shape[0]):
                    # Set field value

                    cmd('XIM F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   xim[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets Y real image height, if provided

    if yri is not None:

        # Update Y real image height for all zoom positions

        if yri.ndim == 1:
            cmd('YRI' + vec_to_str(yri))

        # Update Y real image height for different zoom positions

        elif yri.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(yri.shape[1])
            set_num_z(yri.shape[0])

            # Loop over field

            for f in range(yri.shape[1]):

                # Zoom field position

                cmd('ZOO YRI F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(yri.shape[0]):
                    # Set field value

                    cmd('YRI F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   yri[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets X real image height, if provided

    if xri is not None:

        # Update X real image height for all zoom positions

        if xri.ndim == 1:
            cmd('XRI' + vec_to_str(xri))

        # Update X real image height for different zoom positions

        elif xri.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(xri.shape[1])
            set_num_z(xri.shape[0])

            # Loop over field

            for f in range(xri.shape[1]):

                # Zoom field position

                cmd('ZOO XRI F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(xri.shape[0]):
                    # Set field value

                    cmd('XRI F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(
                        f + 1,
                        z + 1,
                        xri[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets Y object height, if provided

    if yob is not None:

        # Update Y object height for all zoom positions

        if yob.ndim == 1:
            cmd('YOB' + vec_to_str(yob))

        # Update Y object height for different zoom positions

        elif yob.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(yob.shape[1])
            set_num_z(yob.shape[0])

            # Loop over field

            for f in range(yob.shape[1]):

                # Zoom field position

                cmd('ZOO YOB F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(yob.shape[0]):
                    # Set field value

                    cmd('YOB F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   yob[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Sets X object height, if provided

    if xob is not None:

        # Update X object height for all zoom positions

        if xob.ndim == 1:
            cmd('XOB' + vec_to_str(xob))

        # Update X object height for different zoom positions

        elif xob.ndim == 2:

            # Ensure there are enough field and zoom positions

            set_num_f(xob.shape[1])
            set_num_z(xob.shape[0])

            # Loop over field

            for f in range(xob.shape[1]):

                # Zoom field position

                cmd('ZOO XOB F{0:0.0f}'.format(f + 1))

                # Loop over zooms

                for z in range(xob.shape[0]):
                    # Set field value

                    cmd('XOB F{0:0.0f} Z{1:0.0f} {2:0.12f}'.format(f + 1, z + 1,
                                                                   xob[z, f]))

        # > 2 dimensional arrays not supported

        else:
            sys.exit('Invalid input')

    # Set vignetting, if desired

    if set_vig:
        cmd('in cv_macro:setvig')

    return


def create_ana_zoom(ana_zoom):

    """


    Makes in CODE V an anamorphic zoom design based on a first order layout.


    ana_zoom:       anamorphic zoom object

    """

    # Add surfaces

    cmd('INS S1')

    # Set system settings

    set_pupil(fno=ana_zoom.config.fno * ana_zoom.config.ana_rat)
    set_wl(ana_zoom.config.wl, ana_zoom.config.wl_ref)
    set_fld(yan=np.degrees(ana_zoom.config.yan[ana_zoom.ind, :]),
            xan=np.degrees(ana_zoom.config.xan[ana_zoom.ind, :]),
            set_vig=False)

    cmd('LIS')

    # End CODE V session

    stop_cv(cv)

    return

# Initialize CODE V session

cv = init_cv()