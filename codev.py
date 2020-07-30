"""
Author:         David Henry Lippman
File:           codev.py
Date created:   07/21/20
Date modified:  07/21/20

"""

from utilities import find_roc
import numpy as np
from win32com.client import DispatchWithEvents
import pythoncom
import sys


# class DispatchWithEvents:
#
#     def __init__(self, arg1, arg2):
#         self.StartingDirectory = ''
#
#     def StartCodeV(self):
#         return
#
#     def Command(self, arg1):
#         return


class ICVApplicationEvents:

    def OnLicenseError(self, error):
        # This event handler is called when a licensing error is
        # detected in the CODE V application.
        print("License error: %s " % error)

    def OnCodeVError(self, error):
        # This event handler is called when a CODE V error message is issued
        # print("CODE V error: %s " % error)
        return

    def OnCodeVWarning(self, warning):
        # This event handler is called when a CODE V warning message is issued
        # print("CODE V warning: %s " % warning)
        return

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


def cmd(s, out=False):
    
    """
    
    
    Issues CODE V command.
    
    
    s:          command string

    out:        flag for printing output
    
    """

    result = cv.Command(s)
    if out:
        print(s)
        print(result)
    return result


def eva(s, dtype=float):
    
    """
    
    
    Evaluates CODE V database item.
    
    
    s:          evaluation string without database open-close parentheses

    dtype:      data type of output to return
    
    """

    return dtype(cv.EvaluateExpression('(' + s + ')'))


def num_z():

    """


    Gets the number of zoom positions


    """

    return eva('NUM Z', dtype=int)


def num_f():

    """


    Gets the number of field positions


    """

    return eva('NUM F', dtype=int)


def num_s():

    """


    Gets the number of surfaces


    """

    return eva('NUM S', dtype=int)


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


def set_num_z(numz):

    """


    Ensures there are (at least) the specified number of zoom positions. If not,
    adds additional zoom positions.


    numz:           desired number of zoom positions

    """

    while numz > num_z():
        cmd('INS Z{0:0.0f}+1'.format(num_z()))

    return


def set_num_f(numf):

    """


    Ensures there are (at least) the specified number of field positions. If
    not, adds additional field positions.


    numf:           desired number of zoom positions

    """

    while numf > num_f():
        cmd('INS F{0:0.0f}+1'.format(num_f()))

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
            cmd('EPD {0:0.12f}'.format(epd), out=False)

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


def set_tit(sol_type, efx, efy):
    
    """
    
    
    Sets the title for an anamorphic zoom design
    
    
    sol_type:       solution type, e.g. PNPN

    efx:            (m,) array of X effective focal length for m zoom positions

    efy:            (m,) array of Y effective focal length for m zoom positions

    """

    cmd('ZOO TIT')
    for z, (fx, fy) in enumerate(zip(efx, efy)):
        cmd("TIT Z{0:0.0f} '{1}  /  EFX = {2:0.0f}  /  EFY = {3:0.0f}'"
            .format(z + 1, sol_type, fx, fy))


def prv_cat(name, n, wl):
    
    """
    

    Creates a private catalog glass

    
    name:       name of private catalog glass to create

    n:          (m,) array of scalar of refractive index value(s) corresponding
                to the defined wavelengths

    wl:         (m, ) array of wavelengths over which the refractive index
                values are defined [nm]

    """

    # Ensure wavelengths are in an array

    wl = np.array(wl)

    # Ensure refractive index value(s) are in an array

    if np.isscalar(n):
        n *= np.ones_like(wl)

    # Create private catalog

    cmd('PRV')
    cmd('PWL' + vec_to_str(wl))
    cmd("'{0}'".format(name) + vec_to_str(n))
    cmd('END')


def save_seq(filename):

    """


    Save anamorphic zoom design as CODE V sequence file


    filename:       filename to save

    """

    cmd("WRL '{0}'".format(filename))


def create_ana_zoom(ana_zoom, filename=None):

    """


    Makes in CODE V an anamorphic zoom design based on a first order layout.


    ana_zoom:       anamorphic zoom object

    filename:       filename to save model as an SEQ file; if not provided,
                    model is not saved (default)

    """

    # Initialize variables

    num_group_ele = 3

    gla_n_pos = 1.5168  # N-BK7
    gla_name_pos = 'GLA_P'

    gla_n_neg = 1.7552  # N-SF4
    gla_name_neg = 'GLA_N'

    # Reset system

    cmd('in cv_macro:cvnewlens_og')  # note: other users will have cvnewlens.seq

    # Set system settings

    cmd('DIM M')
    set_pupil(fno=ana_zoom.config.fno * ana_zoom.config.ana_rat)
    set_wl(ana_zoom.config.wl, ana_zoom.config.wl_ref)
    set_fld(yan=np.degrees(ana_zoom.config.yan[ana_zoom.ind, :]),
            xan=np.degrees(ana_zoom.config.xan[ana_zoom.ind, :]),
            set_vig=False)
    set_tit(ana_zoom.sol_type, ana_zoom.efx, ana_zoom.efy)

    # Create private catalog glasses

    prv_cat(gla_name_pos, gla_n_pos, ana_zoom.config.wl)
    prv_cat(gla_name_neg, gla_n_neg, ana_zoom.config.wl)

    # Add surfaces

    cmd('INS S2..{0:0.0f}'.format(2 * num_group_ele * ana_zoom.num_group))

    # Set radii of curvature and zoomed air spaces

    # Loop over groups

    for g in range(ana_zoom.num_group):

        # Determine and set surface type

        rad_orient = 'RDY'
        if ana_zoom.group_type[g] == 'X':
            rad_orient = 'RDX'

        # Calculate radii of curvature for all elements in group

        if ana_zoom.group_efl[g] > 0:
            R1, R2 = find_roc(ana_zoom.group_efl[g] * num_group_ele, gla_n_pos)
        else:
            R1, R2 = find_roc(ana_zoom.group_efl[g] * num_group_ele, gla_n_neg)

        # Loop over elements in group

        for e in range(num_group_ele):

            # Calculate front and rear surface numbers

            S1 = g * 2 * num_group_ele + 2 * e + 1
            S2 = g * 2 * num_group_ele + 2 * e + 2

            # Set surface type

            if ana_zoom.group_type[g] == 'X' or ana_zoom.group_type[g] == 'Y':
                cmd('CYL S{0:0.0f}'.format(S1))
                cmd('CYL S{0:0.0f}'.format(S2))

            # Assign radii of curvature

            cmd('{0:3s} S{1:0.0f} {2:0.12f}'.format(rad_orient, S1, R1))
            cmd('{0:3s} S{1:0.0f} {2:0.12f}'.format(rad_orient, S2, R2))

            # Set glass type

            if ana_zoom.group_efl[g] > 0:
                cmd("GLA S{0:0.0f} '{1}'".format(S1, gla_name_pos))
            else:
                cmd("GLA S{0:0.0f} '{1}'".format(S1, gla_name_neg))

            # Add group surface labels

            if e == 0:
                cmd("SLB S{0:0.0f} 'G{1:0.0f}S'".format(S1, g + 1))

            if e == num_group_ele - 1:
                cmd("SLB S{0:0.0f} 'G{1:0.0f}E'".format(S2, g + 1))

        # Zoom airspaces between groups, except for the BFL

        if g < ana_zoom.num_group - 1:

            # Zoom thickness

            cmd("ZOO THI S'G{0:0.0f}E'".format(g + 1))

            # Loop over zoom positions

            for z in range(ana_zoom.num_zoom):

                # Set airspace

                airspace = ana_zoom.group_z[g, z] - ana_zoom.group_z[g + 1, z]

                cmd("THI S'G{0:0.0f}E' Z{1:0.0f} {2:0.12f}".format(g + 1, z + 1,
                                                                   airspace))

    # Set BFL

    cmd("THI S'G{0:0.0f}E' {1:0.12f}".format(ana_zoom.num_group,
                                             ana_zoom.bfl))

    # Add stop surface

    cmd('INS SI-{0:0.0f}'.format(2 * num_group_ele))
    cmd('STO SI-{0:0.0f}'.format(2 * num_group_ele + 1))

    # Save SEQ file, if desired

    if filename is not None:
        save_seq(filename)

    return


def ray_trace():

    """


    Checks for ray trace failures at full field for all defined zoom positions


    """

    # Loop over zoom positions

    for z in range(1, num_z() + 1):

        # Loop over full fields in X, Y and XY only

        for f in range(1, num_f() + 1):

            # Loop over reference rays

            for r in range(1, 6):

                # Trace ray

                cmd('RSI Z{0:0.0f} F{1:0.0f} R{2:0.0f}'.format(z, f, r))

                # Check for ray error

                if eva('RER'):

                    return False

    return True


def opti_ana_zoom():

    """


    Optimizes a first order, ray traceable anamorphic zoom solution. Varies only
    radii of curvature of image defocus.


    """

    # Vary radii of curvature

    # Loop over surfaces

    for s in range(1, num_s()):

        # Spherical surface

        if eva('TYP SUR S{0:0.0f}'.format(s), dtype=str) == 'SPH':

            cmd('CCY S{0:0.0f} 0'.format(s))

        # Cylindrical surface

        else:

            cmd('CCY S{0:0.0f} 0'.format(s))
            cmd('CCX S{0:0.0f} 0'.format(s))

    # Freeze stop surface ROC

    cmd('CCY SS 100')

    # Vary image defocus

    cmd('THC SI 0')

    # Optimize

    opti = 'AUT;'

    # Loop over zoom positions

    for z in range(1, num_z() + 1):

        opti += 'EFY Z{0:0.0f} = {1:0.12f};'\
                .format(z, eva('EFY Z{0:0.0f}'.format(z), dtype=float))
        opti += 'EFX Z{0:0.0f} = {1:0.12f};'\
                .format(z, eva('EFX Z{0:0.0f}'.format(z), dtype=float))

    opti += 'GO'

    cmd(opti)

    return


def avg_spot_size():

    """


    Calculates the average spot size across all zooms and fields for a ray
    traceable anamorphic zoom design


    """

    # Initialize variables

    spo = np.empty((num_f(), num_z()))

    # Calculate RMS spot size for all zooms and fields

    # Initialize data output array for CODE V defined function, SPOTDATA

    cmd('LCL NUM ^arr(10)')

    # Loop over zoom positions

    for z in range(num_z()):

        # Loop over field positions

        for f in range(num_f()):

            # Get spot size from SPOTDATA

            cmd('^dum == SPOTDATA({0:0.0f}, {1:0.0f}, 1, 0, '
                                 '"DEF", 0, 0, ^arr)'.format(z + 1, f + 1))

            spo[f, z] = eva('^arr(1)', dtype=float)

    # Return average spot size

    return spo.mean()


def avg_group_efl(num_group):

    """


    Calculates the average group EFL (absolute value)


    num_group:      number of groups (used for surface labels)

    """

    # Initialize variables

    group_efl = []

    # Loop over groups

    for g in range(1, num_group + 1):

        # Handles both spherical and cylindrical cases

        efx = abs(eva("EFX S'G{0:0.0f}S'..'G{0:0.0f}E'".format(g), dtype=float))
        efy = abs(eva("EFY S'G{0:0.0f}S'..'G{0:0.0f}E'".format(g), dtype=float))

        group_efl.append(min(efx, efy))

    # Return average group EFL (absolute value)

    return np.array(group_efl).mean()


def avg_clear_aper():

    """


    Calculates the average element clear aperture across all surfaces


    """

    # Initialize variables

    ca = np.empty((num_s() - 1, num_z()))

    # Loop over surfaces

    for s in range(num_s() - 1):

        # Loop over zooms

        for z in range(num_z()):

            ca[s, z] = 2 * eva('MAP S{0:0.0f} Z{1:0.0f}'.format(s + 1, z + 1),
                               dtype=float)

    # Return average element clear aperture

    return ca.mean()


def tho():

    # Initialize variables

    tho = np.empty((5, num_z()))

    # Loop over zooms

    for z in range(num_z()):

        # Spherical aberration

        tho[0, z] = eva('SA Z{0:0.0f}'.format(z + 1))

        # Coma

        tho[1, z] = eva('TCO Z{0:0.0f}'.format(z + 1))

        # Astigmatism

        tho[2, z] = eva('TAS Z{0:0.0f}'.format(z + 1))
        tho[3, z] = eva('SAS Z{0:0.0f}'.format(z + 1))

        # Petzval

        tho[4, z] = eva('PTB Z{0:0.0f}'.format(z + 1))

    return tho


# Initialize CODE V session

cv = init_cv()
cmd('REC NO')  # disabling data recording speeds up COM execution
