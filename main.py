"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/20/20

"""

from codev import create_ana_zoom
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def main():

    # Set system values

    config = SystemConfig()

    sols = mc_search_cyl_var(config)

    create_ana_zoom(sols[0])

    """
    
    # STEP 1: system configuration class
    
    # STEP 2: Monte Carlo search
    
    # STEP 3: Anamorphic zoom object
    
    STEP 4: CODE V analysis
    
        4a) Create CODE V design
        
        4b) Check for ray trace failures at all zooms and fields
        
        4c) For successful ray traces, optimize and evluate THO
    
    STEP 5: Demographics
    
        5a) Sankey plot
        
    STEP 6: Spherical variator
    
    """

    return


def mc_search_cyl_var(config, num_trial=1e6):

    """


    Performs a Monte Carlo of cylindrical variator design forms.


    config:         system configuration object

    num_trial:      number of Monte Carlo trials

    """

    # Initialize variables

    num_trial = int(num_trial)

    design_forms = np.array([[ 1,  1,  1,  1],
                             [ 1,  1,  1, -1],
                             [ 1,  1, -1,  1],
                             [ 1,  1, -1, -1],
                             [ 1, -1,  1,  1],
                             [ 1, -1,  1, -1],
                             [ 1, -1, -1,  1],
                             [ 1, -1, -1, -1],
                             [-1,  1,  1,  1],
                             [-1,  1,  1, -1],
                             [-1,  1, -1,  1],
                             [-1,  1, -1, -1],
                             [-1, -1,  1,  1],
                             [-1, -1,  1, -1],
                             [-1, -1, -1,  1],
                             [-1, -1, -1, -1]])

    sols = []

    # Loop over Monte Carlo trials

    for i in range(num_trial):

        # Guess random values for system TTL and BFL

        bfl = rand_rng(config.bfl_rng[0], config.bfl_rng[1], sign=1)
        ttl = rand_rng(config.ttl_rng[0], config.ttl_rng[1], sign=1)
        oal = ttl - bfl

        # Randomly pick a design form

        design_form = design_forms[np.random.randint(0, design_forms.shape[0])]

        # Look for solution(s) in X

        # Guess random group focal lengths including for the stationary groups
        # and x-oriented moving groups

        f1 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[0])
        f2_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[1])
        f3_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[2])
        f4 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[3])

        # Check for valid four group zoom solution(s)

        sols_x = check_four_group(config.efx, oal, bfl, f1, f2_x, f3_x, f4,
                                  min_air=10)

        # Loop over solutions in X, if any were found

        for sol_x in sols_x:

            # Look for solution(s) in Y

            # Guess random group focal lengths including for the y-oriented
            # moving groups only, not the stationary first and fourth groups
            # which are the same as in the X solution

            f2_y = rand_rng(config.efl_group_rng[0],
                            config.efl_group_rng[1], sign=design_form[1])
            f3_y = rand_rng(config.efl_group_rng[0],
                            config.efl_group_rng[1], sign=design_form[2])

            # Check for valid four group zoom solution(s) taking into
            # account solution found in X

            sols_y = check_four_group(config.efy, oal, bfl, f1, f2_y, f3_y,
                                      f4, min_air=10, sol_check=sol_x)

            # Loop over solutions in X and Y, if any were found

            for sol_y in sols_y:

                # Combine X and Y solutions for group focal lengths, zoom
                # motions, and group types

                group_efl, \
                group_z, \
                group_type = combine_sols(sol_x, sol_y, f1, f2_x, f2_y, f3_x,
                                          f3_y, f4)

                # Create anamorphic zoom solution and add to solution
                # list

                ana_zoom = AnamorphicZoom(config, group_efl, group_type,
                                          group_z)
                print(ana_zoom)

                # Plot solution

                ana_zoom.plot_zoom()

                # Add to solution list

                sols.append(ana_zoom)

                return sols

    return sols


def rand_rng(min_val, max_val, sign=1):

    """


    Returns a random value within a range with the desired sign


    min_val:        minimum value of range (absolute value)

    max_val:        maximum value of range (absolute value)

    sign:           the sign of the random value;

                    1 for positive
                   -1 for negative
                    0 for random

    """

    if sign == 0 and np.random.random() > 0.5:
        sign = 1
    elif sign == 0:
        sign = -1

    return sign * np.random.uniform(np.abs(min_val), np.abs(max_val))


def check_four_group(efl, oal, bfl, f1, f2, f3, f4, min_air=0., sol_check=None):

    """


    Identifies whether there is a valid solution for a four group zoom design.
    Checks both positive and negative roots.

    References:

    [1] A. J. Yee, D. J. L. Williams, G. A. Gandara-Montano, P. McCarthy,
        J. Bentley, and D. T. Moore, “New tools for finding first-order zoom
        lens solutions and the analysis of zoom lenses during the design
        process,” San Diego, California, United States, Sep. 2015, p. 958006,
        doi: 10.1117/12.2186780.


    efl:            effective focal lengths of zoom solution [mm]

    oal:            overall length of zoom solution [mm]

    bfl:            back focal length of zoom solution [mm]

    f1:             group 1 focal length [mm]

    f2:             group 2 focal length [mm]

    f3:             group 3 focal length [mm]

    f4:             group 4 focal length [mm]

    min_air:        minimum air space between groups; default is 0

    sol_check:      previously found solution to check with for crashes

    """

    # Initialize data structure

    sols = []

    # Calculate quadratic equation coefficients

    s4 = -f4 * bfl / (bfl - f4)
    L = oal - f1 + s4
    M = -efl * (oal - L - f1) / (bfl * f1)
    c = L * (f2 + f3) + ((M - 1) ** 2) * f2 * f3 / M

    # Check if roots are real for all zoom positions (all values of c)

    if L ** 2 - 4 * c.min() >= 0 and L ** 2 - 4 * c.max() >= 0:

        # Check positive root

        t2 = (L + np.sqrt(L ** 2 - 4 * c)) / 2

        # Calculate zoom motions

        z = calc_zoom_motion(M, L, t2, oal, bfl, f2, s4, efl)

        # Check if there is a valid solution (no crashes)

        if np.diff(z, axis=0).max() < -min_air:

            # Check other solution, if provided

            if sol_check is not None:

                # Look at zoom motions of moving groups and sort by arbitrary
                # (first) zoom position

                z_moving = np.concatenate((z[1: -1, :], sol_check[1: -1, :]))
                z_moving = z_moving[z_moving[:, 0].argsort()[::-1]]

                # Check if there is a valid solution (no crashes)

                if np.diff(z_moving, axis=0).max() < -min_air:

                    # Add solution to array

                    sols.append(z)

            # Add solution to array

            else:
                sols.append(z)

        # Check negative root

        t2 = (L - np.sqrt(L ** 2 - 4 * c)) / 2

        if t2.min() > 0:

            # Calculate zoom motions

            z = calc_zoom_motion(M, L, t2, oal, bfl, f2, s4, efl)

            # Check if there is a valid solution (no crashes)

            if np.diff(z, axis=0).max() < -min_air:

                # Check other solution, if provided

                if sol_check is not None:

                    # Look at zoom motions of moving groups and sort by arbitrary
                    # (first) zoom position

                    z_moving = np.concatenate((z[1: -1, :],
                                               sol_check[1: -1, :]))
                    z_moving = z_moving[z_moving[:, 0].argsort()[::-1]]

                    # Check if there is a valid solution (no crashes)

                    if np.diff(z_moving, axis=0).max() < -min_air:

                        # Add solution to array

                        sols.append(z)

                # Add solution to array

                else:
                    sols.append(z)



    return sols


def calc_zoom_motion(M, L, t2, oal, bfl, f2, s4, efl):

    """


    Finds the zoom motions for a four group zoom.

    References:

    [1] A. J. Yee, D. J. L. Williams, G. A. Gandara-Montano, P. McCarthy,
        J. Bentley, and D. T. Moore, “New tools for finding first-order zoom
        lens solutions and the analysis of zoom lenses during the design
        process,” San Diego, California, United States, Sep. 2015, p. 958006,
        doi: 10.1117/12.2186780.


    M:              M quadratic coefficient

    L:              L quadratic coefficient

    t2:             air space between groups 2 and 3 [mm]

    oal:            overall length of zoom solution [mm]

    bfl:            back focal length of zoom solution [mm]

    f2:             group 2 focal length [mm]

    s4:             front focal distance of group 4 [mm]

    efl:            effective focal lengths of zoom solution [mm]

    """

    # Calculate air spaces

    s2 = ((M - 1) * t2 + L) / ((M - 1) - M * t2 / f2)
    t1 = oal - L + s4 - s2
    t3 = oal - t1 - t2

    # Calculate group zoom motions

    z = np.ones((4, efl.size)) * bfl
    z[2, :] = t3 + z[3, :]
    z[1, :] = t2 + z[2, :]
    z[0, :] = t1 + z[1, :]

    return z


def combine_sols(sol_x, sol_y, f1, f2_x, f2_y, f3_x, f3_y, f4):

    # Initialize data structures

    group_efl = np.empty(6)
    group_z = np.empty((6, sol_x.shape[1]))
    group_type = np.empty(6, dtype='<U2')

    # Store values for stationary first and last groups

    group_efl[0] = f1
    group_efl[-1] = f4

    group_z[0, :] = sol_x[0, :]
    group_z[-1, :] = sol_x[-1, :]

    group_type[0] = 'XY'
    group_type[-1] = 'XY'

    # Determine variator order in X and Y

    if sol_x[1, 0] > sol_y[1, 0]:
        group_efl[1] = f2_x
        group_efl[2] = f2_y
        group_z[1, :] = sol_x[1, :]
        group_z[2, :] = sol_y[1, :]
        group_type[1] = 'X'
        group_type[2] = 'Y'
    else:
        group_efl[1] = f2_y
        group_efl[2] = f2_x
        group_z[1, :] = sol_y[1, :]
        group_z[2, :] = sol_x[1, :]
        group_type[1] = 'Y'
        group_type[2] = 'X'

    # Determine compensator order in X and Y

    if sol_x[2, 0] > sol_y[2, 0]:
        group_efl[3] = f3_x
        group_efl[4] = f3_y
        group_z[3, :] = sol_x[2, :]
        group_z[4, :] = sol_y[2, :]
        group_type[3] = 'X'
        group_type[4] = 'Y'
    else:
        group_efl[3] = f3_y
        group_efl[4] = f3_x
        group_z[3, :] = sol_y[2, :]
        group_z[4, :] = sol_x[2, :]
        group_type[3] = 'Y'
        group_type[4] = 'X'

    return group_efl, group_z, group_type


class SystemConfig:

    """


    Anamorphic zoom system setup.


    """

    def __init__(self, ana_rat=2., efl_sys_rng=np.array([28., 76.]),
                 bfl_rng=np.array([35., 65.]), ttl_rng=np.array([240., 365.]),
                 efl_group_rng=np.array([20, 500]), num_zoom=45, fno=4.,
                 img_dim=np.array([22.31, 18.67]),
                 wl=np.array([656.3, 587.6, 486.1]), wl_ref=1):
        
        """


        Initialize anamorphic zoom system object.

        
        ana_rat:            anamorphic ratio

        efl_sys_rng:        (2,) array of system effective focal length range
                            [mm]; defined in horizontal a.k.a. axis

        bfl_rng:            (2,) array of back focal length range [mm]

        ttl_rng:            (2,) array of total track length range [mm]; defined
                            from front surface vertex to image

        efl_group_rng:      (2,) array of group effective focal length range
                            [mm]; really the absolute value

        num_zoom:           number of zoom positions to evaluate in the design

        fno:                f/number

        img_dim:            (2,) array of image dimensions [mm]

        wl:                 (m,) array of wavelength(s) [nm]

        wl_ref:             reference wavelength index (zero-indexed unlike
                            CODE V which is one-indexed)

        """

        # Initialize variables

        # Monte Carlo search

        self.ana_rat = ana_rat
        self.efl_sys_rng = efl_sys_rng
        self.bfl_rng = bfl_rng
        self.ttl_rng = ttl_rng
        self.efl_group_rng = efl_group_rng

        # System

        self.num_zoom = num_zoom
        self.fno = fno
        self.img_dim = img_dim
        self.wl = wl
        self.wl_ref = wl_ref
        self.num_fld = 13

        # Calculate EFL in Y and X for all zoom positions

        self.efx = np.linspace(self.efl_sys_rng[0], self.efl_sys_rng[1],
                               self.num_zoom)
        self.efy = self.efx * 2

        # Set field of view

        self.set_fov()

    def __repr__(self):

        # TODO

        repr_str = ''

        return repr_str

    def set_fov(self, scale_fact=1.):

        """


        Sets the system angular field-of-view.


        scale_fact:     scale factor for image space field of view; default is 1

        """

        # Form proportional image space FOV vectors for unity HFOV

        xim_prop = np.array([0.,
                             0.,
                             0.,
                             0.,
                             0.,
                             0.5,
                             np.sqrt(2) / 2,
                             np.sqrt(3) / 2,
                             1.,
                             0.5,
                             np.sqrt(2) / 2,
                             np.sqrt(3) / 2,
                             1.])

        yim_prop = np.array([0.,
                             0.5,
                             np.sqrt(2) / 2,
                             np.sqrt(3) / 2,
                             1.,
                             0.,
                             0.,
                             0.,
                             0.,
                             0.5,
                             np.sqrt(2) / 2,
                             np.sqrt(3) / 2,
                             1.])

        # Calculate image space FOV

        xim = xim_prop * self.img_dim[0] / 2
        yim = yim_prop * self.img_dim[1] / 2

        # Calculate angular FOV

        self.xan = np.arctan(xim / np.reshape(self.efx, (self.num_zoom, 1)))
        self.yan = np.arctan(yim / np.reshape(self.efy, (self.num_zoom, 1)))

        # Scale by scale factor

        self.xan *= scale_fact
        self.yan *= scale_fact

        return

    def plot_fov(self):

        """


        Plots the field-of-view.


        """

        # Loop over zoom positions

        for z in [0, -1]:

            # Initialize figure

            fig, ax = plt.subplots()

            # Plot field points

            ax.scatter(np.degrees(self.xan[:, z]),
                       np.degrees(self.yan[:, z]))

            # Plot settings

            ax.set(title='EFX = {0:0.2f} mm\nEFY = {1:0.2f} mm'
                         .format(self.efx[z], self.efy[z]),
                   xlabel='X field angle [deg]', ylabel='Y field angle [deg]')
            ax.set_aspect('equal')


class AnamorphicZoom:

    """


    Anamorphic zoom design object


    """

    def __init__(self, config, group_efl, group_type, group_z, num_zoom=5):

        """


        Initialize anamorphic zoom design object


        config:         system configuration object

        group_efl:      (m,) array of group effective focal lengths [mm]

        group_type:     (m,) array of group types [str];
                        'XY' = spherical
                        'X'  = cylindrical oriented in X
                        'Y'  = cylindrical oriented in Y

        group_z:        (m, n) array of zoom motions for the n evaluated zoom
                        positions [mm]

        num_zoom:       number of zoom positions to consider for the design

        """

        # Initialize variables

        self.config = config
        self.group_efl = group_efl
        self.group_type = group_type
        self.group_z = group_z
        self.num_zoom = num_zoom

        # Calculate basic values

        self.num_group = self.group_efl.size
        self.ttl = self.group_z[0, 0]
        self.bfl = self.group_z[-1, 0]
        self.oal = self.ttl - self.bfl

        # Calculate even spacing indices

        self.ind = np.linspace(0, config.num_zoom - 1, 5).astype(int)

        # Calculate zoom position focal lengths

        self.efx = config.efx[self.ind]
        self.efy = config.efy[self.ind]

        # Classify solution type

        self.classify()

    def __repr__(self):

        repr_str = '\n' + 10 * '~' + ' ANAMORPHIC ZOOM SOLUTION ' + 10 * '~' + \
                   '\n\n'

        repr_str += 'Variator type:\t\t{0}\n'.format(self.vari_str.capitalize())
        repr_str += 'Solution type:\t\t{0}\n'.format(self.sol_str)
        repr_str += 'Cyl. orient.:\t\t{0}\n\n'.format(self.orient_str)

        repr_str += 'Groups:\n\n'

        template = '{0:>4.0f}{1:>8s}{2:>8s}{3:>14.4f} mm\n'
        for i in range(self.num_group):

            group = np.ceil(i / 2).astype(int) + 1
            surf_type = 'CYL'
            if self.group_type[i] == 'XY':
                surf_type = 'SPH'

            surf_orient = '-'
            if surf_type == 'CYL':
                surf_orient = self.group_type[i]

            repr_str += template.format(group, surf_type, surf_orient,
                                        self.group_efl[i])

        repr_str += '\n' + (2 * 10 + 26) * '~' + '\n\n'

        return repr_str

    def classify(self):

        # Classify power

        self.sol_str = ''

        for i in [0, 1, 3, -1]:
            if self.group_efl[i] > 0:
                self.sol_str += 'P'
            else:
                self.sol_str += 'N'

        # Classify cylinder orientation

        self.orient_str = ''
        for char in self.group_type[self.group_type != 'XY']:
            self.orient_str += char

        # Classify variator type

        if self.num_group == 5:
            self.vari_str = 'spherical'
        elif self.num_group == 6:
            self.vari_str = 'cylindrical'

    def plot_zoom(self):

        # Initialize plot variables

        clrs = np.empty_like(self.group_type)
        clrs[self.group_type == 'XY'] = 'b'
        clrs[self.group_type == 'X'] = 'g'
        clrs[self.group_type == 'Y'] = 'g'

        lstyles = np.empty_like(self.group_type)
        lstyles[self.group_type == 'XY'] = '-'
        lstyles[self.group_type == 'X'] = '-'
        lstyles[self.group_type == 'Y'] = '--'

        # Initialize figure

        fig, ax = plt.subplots(figsize=(12, 6))
        plt.subplots_adjust(left=0.25, right=0.65)

        # Create second y-axis

        ax_par = ax.twinx()
        ax_par.spines["left"].set_position(("axes", -0.2))
        ax_par.set_frame_on(True)
        ax_par.patch.set_visible(False)
        for sp in ax_par.spines.values():
            sp.set_visible(False)
        ax_par.spines["left"].set_visible(True)
        ax_par.yaxis.set_label_position('left')
        ax_par.yaxis.set_ticks_position('left')

        # Plot zoom

        for i in range(self.num_group):

            # Form legend string

            group_ind = np.ceil(i / 2).astype(int) + 1

            type_str = self.group_type[i]
            if type_str == 'XY':
                type_str = ''

            subscript_str = '{0}'.format(group_ind) + type_str
            leg_str = r'$f_{{{0:2s}}}$ = {1:0.2f} mm'.format(subscript_str,
                                                             self.group_efl[i])

            # Plot zoom motion

            ax.plot(self.group_z[i, :], self.config.efx, color=clrs[i],
                    linestyle=lstyles[i], label=leg_str)

            ax_par.plot(self.group_z[i, :], self.config.efy, color='none')

        # Set plot settings

        ax.set(title=self.sol_str + '  /  ' + self.orient_str,
               xlabel="z [mm]",
               ylabel="EFX [mm]")
        ax_par.set(ylabel="EFY [mm]")

        ax.set_xlim(self.ttl * 1.1, 0)

        ax.set_yticks(self.efx)
        ax.set_yticklabels(self.efx)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax_par.set_yticks(self.efy)
        ax_par.set_yticklabels(self.efy)
        ax_par.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

        return


if __name__ == '__main__':
    main()
    plt.show()

