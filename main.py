"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/20/20

"""

import numpy as np
import matplotlib.pyplot as plt


def main():

    # Set system values

    config = SystemConfig()

    mc_search_cyl_var(config)

    return


def mc_search_cyl_var(config, num_trial=1e6):

    """


    Performs a Monte Carlo of cylindrical variator design forms.


    config:         system configuration object

    num_trial:      number of Monte Carlo trials

    """

    # design_form = np.array([[ 1,  1,  1,  1],
    #                         [ 1,  1,  1, -1],
    #                         [ 1,  1, -1,  1],
    #                         [ 1,  1, -1, -1],
    #                         [ 1, -1,  1,  1],
    #                         [ 1, -1,  1, -1],
    #                         [ 1, -1, -1,  1],
    #                         [ 1, -1, -1, -1],
    #                         [-1,  1,  1,  1],
    #                         [-1,  1,  1, -1],
    #                         [-1,  1, -1,  1],
    #                         [-1,  1, -1, -1],
    #                         [-1, -1,  1,  1],
    #                         [-1, -1,  1, -1],
    #                         [-1, -1, -1,  1],
    #                         [-1, -1, -1, -1]])

    # Initialize variables

    design_form = np.array([1, -1, 1, 1])

    # TEMP: Loop until a valid solution is found

    flag = True

    while flag:

        # Guess random values for system TTL and BFL

        bfl = rand_rng(config.bfl_rng[0], config.bfl_rng[1], sign=1)
        ttl = rand_rng(config.ttl_rng[0], config.ttl_rng[1], sign=1)
        oal = ttl - bfl

        # STEP 1: Look for solution(s) in X

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

        # Solution(s) found in X

        if len(sols_x) > 0:

            # Loop over solutions in X

            for sol_x in sols_x:

                # STEP 2: Look for solution(s) in Y

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

                # Solution(s) found in Y

                if len(sols_y) > 0:

                    # Loop over solutions in X and Y

                    for sol_y in sols_y:

                        # Plot solution

                        fig, ax = plt.subplots()

                        ax.plot(config.efx, sol_x.transpose(),
                                color='blue')
                        ax.plot(config.efx, sol_y[1: -1, :].transpose(),
                                color='green')

                        flag = False

    return


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


class SystemConfig:

    """


    Anamorphic zoom system setup.


    """

    def __init__(self, ana_rat=2., efl_sys_rng=np.array([28., 76.]),
                 bfl_rng=np.array([35., 65.]), ttl_rng=np.array([240., 365.]),
                 efl_group_rng=np.array([20, 500]), num_zoom=5, fno=4.,
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

        num_zoom:           number of zoom positions

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

        self.xan = np.arctan(np.reshape(xim, (self.num_fld, 1)) / self.efx)
        self.yan = np.arctan(np.reshape(yim, (self.num_fld, 1)) / self.efy)

        # Scale by scale factor

        self.xan *= scale_fact
        self.yan *= scale_fact

        # # TEMP: plot
        #
        # fig, axs = plt.subplots(1, self.num_zoom)
        #
        # for z in range(self.num_zoom):
        #     axs[z].scatter(np.degrees(self.xan[:, z]),
        #                    np.degrees(self.yan[:, z]))
        #     axs[z].set(xlabel='X field angle [deg]',
        #                ylabel='Y field angle [deg]',
        #                title='EFX = {0:0.2f} mm\nEFY = {1:0.2f} mm'
        #                      .format(self.efx[z], self.efy[z]))
        #     axs[z].set_aspect('equal')

        return


class AnamorphicZoom:

    def __init__(self):

        return

if __name__ == '__main__':
    main()
    plt.show()

