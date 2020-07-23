"""
Author:         David Henry Lippman
File:           monte_carlo.py
Date created:   07/22/20
Date modified:  07/22/20

"""

from designs import AnamorphicZoom, Solutions
from codev import opti_ana_zoom, avg_spot_size
from utilities import rand_rng, get_time_str, folder_exist, save_obj
import numpy as np


def mc_search_cyl_var(config, num_trial=1e6):

    """


    Performs a Monte Carlo of cylindrical variator design forms.


    config:         system configuration object

    num_trial:      number of Monte Carlo trials

    """

    # Initialize variables

    num_trial = int(num_trial)
    min_air = 0

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

    sols = Solutions(config, num_trial)

    progress = (num_trial * np.linspace(0.1, 1, 10)).astype(int)

    # Create folder to store solutions ins

    start_time = get_time_str()
    path = 'C:\\CVUSER\\Anamorphic Zoom Solutions\\Cylindrical Variator\\' \
           + start_time + '\\'

    folder_exist(path)

    # Loop over Monte Carlo trials

    print('\nPerforming a Monte Carlo search with {0:0.1e} trials\n'
          .format(num_trial))

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
                                  min_air=min_air)

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
                                      f4, min_air=min_air, sol_check=sol_x)

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
                # print(ana_zoom)

                # Create model in CODE V

                ana_zoom.make_codev()

                # Create model and check ray trace in CODE V

                ana_zoom.check_ray_trace()

                # Successful ray trace

                if not ana_zoom.ray_fail:

                    # Optimize

                    opti_ana_zoom()

                    # Calculate average spot size

                    ana_zoom.avg_spot_size = avg_spot_size()

                    # Save CODE V model

                    ana_zoom.save_seq(sol_num=i, path=path)

                # Add to general solution list

                sols.add_sol(ana_zoom)

        # Check status

        if i + 1 in progress:
            perc = 100 * (np.asscalar(np.argwhere(i + 1 == progress)) + 1) / \
                   progress.size
            print('{0:0.0f}% complete'.format(perc))
            print(sols)

    # Save solutions to file

    save_obj(sols, filename=start_time)

    return sols


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

    group_efl = np.array([f1, f2_x, f3_x, f2_y, f3_y, f4])
    group_z = np.concatenate((sol_x[0: -1, :], sol_y[1:, :]))
    group_type = np.array(['XY', 'X', 'X', 'Y', 'Y', 'XY'])

    # Determine variator and compensator order in X and Y by sorting descending
    # z positions for an arbitrary (first) zoom position

    order = group_z[1: -1, 0].argsort()[::-1] + 1

    # Re-order data structures

    group_efl[1: -1] = group_efl[order]
    group_z[1: -1, :] = group_z[order, :]
    group_type[1: -1] = group_type[order]

    return group_efl, group_z, group_type

