"""
Author:         David Henry Lippman
File:           monte_carlo.py
Date created:   07/22/20
Date modified:  07/28/20

"""

from designs import AnamorphicZoom, Solutions
from utilities import rand_rng, get_time_str, folder_exist, format_time_str, \
    save_obj
import numpy as np
from time import time
import sys
import matplotlib.pyplot as plt


def mc_search_cyl_var(config, num_trial=1e6, same_xy=True):

    """


    Performs a Monte Carlo of cylindrical variator design forms.


    config:         system configuration object

    num_trial:      number of Monte Carlo trials

    same_xy:        flag for whether the design form should be the same in X and
                    Y

    """

    # Initialize variables

    num_trial = int(num_trial)
    min_air = 0

    # Create solutions object

    sols = Solutions(config, num_trial, "CYL", same_xy)

    # Create folder to store CODE V solutions in

    time_str = get_time_str()
    path = 'C:\\CVUSER\\Anamorphic Zoom Solutions\\Cylindrical Variator\\' \
           + time_str + '\\'
    folder_exist(path)

    # Set up progress bar and timer

    num_prog = 100
    progress = (num_trial * np.linspace(0.1, 1, num_prog)).astype(int)
    start_time = time()

    # Loop over Monte Carlo trials

    print('\nPerforming a Monte Carlo search with {0:0.1e} trials\n'
          .format(num_trial))

    for i in range(num_trial):

        # Randomly pick values for system TTL and BFL

        ttl = rand_rng(config.ttl_rng[0], config.ttl_rng[1], sign=1)
        bfl = rand_rng(config.bfl_rng[0], config.bfl_rng[1], sign=1)
        oal = ttl - bfl

        # Randomly pick a design form based on whether X and Y have the same
        # solution type

        design_form = 2 * np.random.randint(0, 2, size=(6,)) - 1
        if same_xy:
            design_form[2] = design_form[1]
            design_form[4] = design_form[3]

        # Randomly pick group focal lengths

        f1 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[0])

        f2_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[1])
        f2_y = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[2])

        f3_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[3])
        f3_y = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[4])

        f4 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[5])

        # Check for valid six group zoom solution(s)

        sols_fo = check_cyl_var_sol(config.efx, config.efy, oal, bfl, f1, f2_x,
                                    f2_y, f3_x, f3_y, f4, min_air)

        # Loop over found solution(s), if any

        for sol in sols_fo:

            """

            Successful first order solution

            """

            # Create anamorphic zoom design

            ana_zoom = AnamorphicZoom(config, sol, num_zoom=5,
                                      sol_num=sols.num_sol)

            # Create model in CODE V

            ana_zoom.make_codev()

            # Check ray trace in CODE V at all zooms and fields

            ana_zoom.check_ray_trace()

            # Successful ray trace at all zooms and fields

            if ana_zoom.ray_trace:

                """

                Successful ray traceable solution

                """

                # Optimize first order design

                ana_zoom.optimize()

                # Analyze optimized CODE V design for performance,
                # aberrations, and packaging

                ana_zoom.analyze()

                # Check that optimization converged

                if ana_zoom.avg_spot_size > 5:  # mm

                    # If not, ray tracing failed (some rays trace
                    # successfully but are behaving erratically

                    ana_zoom.ray_trace = False

                else:

                    # Write design to CODE V sequence file

                    ana_zoom.save_seq(path=path)

            # Add design to solutions

            sols.add_sol(ana_zoom)

        # Check status

        if i + 1 in progress:

            perc = (np.asscalar(np.argwhere(i + 1 == progress)) + 1) / num_prog

            cur_time = time()
            pass_time = cur_time - start_time
            full_time = pass_time / perc
            remain_time = full_time - pass_time

            print('{0:0.0f}% complete\n'.format(perc * 100))
            print('Approximately {0} remaining'
                  .format(format_time_str(remain_time)))
            print(sols)

    # Close out solution timer

    sols.stopwatch(end=True)

    # Save solutions to file if any are found

    if sols.num_sol:
        save_obj(sols, filename=time_str)

    # Delete folder if no ray traceable solutions were found

    if not sols.num_sol_rt:
        folder_exist(path, erase=True)

    return sols


def check_cyl_var_sol(efx, efy, oal, bfl, f1, f2_x, f2_y, f3_x, f3_y, f4,
                      min_air=0.):

    """


    Identifies whether there is a valid solution for a six group cylindrical
    variator anamorphic zoom design.

    References:

    [1]  A. J. Yee, D. J. L. Williams, G. A. Gandara-Montano, P. McCarthy,
         J. Bentley, and D. T. Moore, “New tools for finding first-order zoom
         lens solutions and the analysis of zoom lenses during the design
         process,” San Diego, California, United States, Sep. 2015, p. 958006,
         doi: 10.1117/12.2186780.



    efx:            effective focal lengths in X of zoom system [mm]

    efy:            effective focal lengths in Y of zoom system [mm]

    oal:            overall length of zoom system [mm]

    bfl:            back focal length of zoom system [mm]

    f1:             group 1 focal length [mm]

    f2_x:           group 2 focal length in X [mm]

    f2_y:           group 2 focal length in Y [mm]

    f3_x:           group 3 focal length in X [mm]

    f3_y:           group 3 focal length in Y [mm]

    f4:             group 4 focal length [mm]

    min_air:        minimum air space between groups; default is 0

    """

    # Initialize data structure

    sols = []

    # Calculate quadratic equation coefficients

    s4 = -f4 * bfl / (bfl - f4)
    L = oal - f1 + s4

    M_x = -efx * (oal - L - f1) / (bfl * f1)
    M_y = -efy * (oal - L - f1) / (bfl * f1)

    c_x = L * (f2_x + f3_x) + ((M_x - 1) ** 2) * f2_x * f3_x / M_x
    c_y = L * (f2_y + f3_y) + ((M_y - 1) ** 2) * f2_y * f3_y / M_y

    # Check if roots are real for all zoom positions (all values of c) for both
    # the X and Y slices

    if L ** 2 - 4 * c_x.min() >= 0 and L ** 2 - 4 * c_x.max() >= 0 and \
       L ** 2 - 4 * c_y.min() >= 0 and L ** 2 - 4 * c_y.max() >= 0:

        # Calculates both roots for t2 in both X and Y

        t2_x_roots = np.array([(L + np.sqrt(L ** 2 - 4 * c_x)) / 2,
                               (L - np.sqrt(L ** 2 - 4 * c_x)) / 2])
        t2_y_roots = np.array([(L + np.sqrt(L ** 2 - 4 * c_y)) / 2,
                               (L - np.sqrt(L ** 2 - 4 * c_y)) / 2])

        # Remove negative valued result, if any

        t2_x_roots = t2_x_roots[np.min(t2_x_roots, axis=1) > 0]
        t2_y_roots = t2_y_roots[np.min(t2_y_roots, axis=1) > 0]

        # Loop over all combinations of positive t2 values

        for t2_x in t2_x_roots:

            for t2_y in t2_y_roots:

                # Calculate zoom motions in X and Y

                group_z_x = calc_zoom_motion(oal, bfl, M_x, t2_x, f2_x, L, s4)
                group_z_y = calc_zoom_motion(oal, bfl, M_y, t2_y, f2_y, L, s4)

                # Combine motions in X and Y

                group_z, \
                group_efl, \
                group_type = combine_cyl_sols(group_z_x, group_z_y, f1, f2_x, f2_y,
                                              f3_x, f3_y, f4)

                # Check for a valid solution (no group crashes)

                if np.diff(group_z[::-1], axis=0).min() > min_air:

                    # Add solution to array

                    sols.append([group_z, group_efl, group_type])

    return sols


def calc_zoom_motion(oal, bfl, M, t2, f2, L, s4):

    """


    Finds the zoom motions for a four group zoom.

    References:

    [1] A. J. Yee, D. J. L. Williams, G. A. Gandara-Montano, P. McCarthy,
        J. Bentley, and D. T. Moore, “New tools for finding first-order zoom
        lens solutions and the analysis of zoom lenses during the design
        process,” San Diego, California, United States, Sep. 2015, p. 958006,
        doi: 10.1117/12.2186780.


    oal:            overall length of zoom solution [mm]

    bfl:            back focal length of zoom solution [mm]

    M:              M quadratic coefficient

    L:              L quadratic coefficient

    t2:             air space between groups 2 and 3 [mm]

    f2:             group 2 focal length [mm]

    s4:             front focal distance of group 4 [mm]

    """

    # Calculate air spaces

    s2 = ((M - 1) * t2 + L) / ((M - 1) - M * t2 / f2)
    t1 = oal - L + s4 - s2
    t3 = oal - t1 - t2

    # Calculate group zoom motions

    z = np.ones((4, M.size)) * bfl
    z[2, :] = t3 + z[3, :]
    z[1, :] = t2 + z[2, :]
    z[0, :] = t1 + z[1, :]

    return z


def combine_cyl_sols(sol_x, sol_y, f1, f2_x, f2_y, f3_x, f3_y, f4):

    """


    Combines solutions in X and Y into a single XY design


    sol_x:      (m, n) array of X solution zoom motion with m groups and n zoom
                positions [mm]

    sol_y:      (m, n) array of Y solution zoom motion with m groups and n zoom
                positions [mm]

    f1:         focal length of group 1 in X and Y [mm]

    f2_x:       focal length of group 2 in X [mm]

    f2_y:       focal length of group 2 in Y [mm]

    f3_x:       focal length of group 3 in X [mm]

    f3_y:       focal length of group 3 in Y [mm]

    f4:         focal length of group 4 in X and Y [mm]

    """

    # Initialize data structures

    group_z = np.concatenate((sol_x[0: -1, :], sol_y[1:, :]))
    group_efl = np.array([f1, f2_x, f3_x, f2_y, f3_y, f4])
    group_type = np.array(['XY', 'X', 'X', 'Y', 'Y', 'XY'])

    # Determine variator and compensator order in X and Y by sorting descending
    # z positions based on an arbitrary (first) zoom position

    order = group_z[1: -1, 0].argsort()[::-1] + 1

    # Re-order data structures

    group_z[1: -1, :] = group_z[order, :]
    group_efl[1: -1] = group_efl[order]
    group_type[1: -1] = group_type[order]

    return group_z, group_efl, group_type


def mc_search_sph_var(config, num_trial=1e6):

    """


    Performs a Monte Carlo of spherical variator design forms.


    config:         system configuration object

    num_trial:      number of Monte Carlo trials

    """

    # Initialize variables

    num_trial = int(num_trial)
    min_air = 0
    same_xy = True  # necessarily since ana_rat = f3_x / f3_y

    # Create solutions object

    sols = Solutions(config, num_trial, "SPH", same_xy)

    # Create folder to store CODE V solutions in

    time_str = get_time_str()
    path = 'C:\\CVUSER\\Anamorphic Zoom Solutions\\Cylindrical Variator\\' \
           + time_str + '\\'
    folder_exist(path)

    # Set up progress bar and timer

    num_prog = 100
    progress = (num_trial * np.linspace(0.1, 1, num_prog)).astype(int)
    start_time = time()

    # Loop over Monte Carlo trials

    print('\nPerforming a Monte Carlo search with {0:0.1e} trials\n'
          .format(num_trial))

    for i in range(num_trial):

        # Randomly pick values for system TTL and BFL

        ttl = rand_rng(config.ttl_rng[0], config.ttl_rng[1], sign=1)
        bfl = rand_rng(config.bfl_rng[0], config.bfl_rng[1], sign=1)
        oal = ttl - bfl

        # Randomly pick a design form
        #
        # Relay group must be positive since f4 = bfl > 0

        design_form = 2 * np.random.randint(0, 2, size=(4,)) - 1
        design_form[-1] = 1

        # Randomly pick group focal lengths

        f1 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[0])

        f2 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[1])

        f3_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[2])
        f3_y = f3_x / config.ana_rat

        f4 = bfl

        # Check for valid five group zoom solution(s)

        sols_fo = check_sph_var_sol(config.efx, oal, bfl, f1, f2, f3_x, f3_y,
                                    f4, min_air)

        # Loop over found solution(s), if any

        for sol in sols_fo:

            """

            Successful first order solution

            """

            # Create anamorphic zoom design

            ana_zoom = AnamorphicZoom(config, sol, num_zoom=5,
                                      sol_num=sols.num_sol)

            # Create model in CODE V

            ana_zoom.make_codev()

            # Check ray trace in CODE V at all zooms and fields

            ana_zoom.check_ray_trace()

            # Successful ray trace at all zooms and fields

            if ana_zoom.ray_trace:

                """

                Successful ray traceable solution

                """

                # Optimize first order design

                ana_zoom.optimize()

                # Analyze optimized CODE V design for performance,
                # aberrations, and packaging

                ana_zoom.analyze()

                # Check that optimization converged

                if ana_zoom.avg_spot_size > 5:  # mm

                    # If not, ray tracing failed (some rays trace
                    # successfully but are behaving erratically

                    ana_zoom.ray_trace = False

                else:

                    # Write design to CODE V sequence file

                    ana_zoom.save_seq(path=path)

            # Add design to solutions

            sols.add_sol(ana_zoom)

        # Check status

        if i + 1 in progress:

            perc = (np.asscalar(np.argwhere(i + 1 == progress)) + 1) / num_prog

            cur_time = time()
            pass_time = cur_time - start_time
            full_time = pass_time / perc
            remain_time = full_time - pass_time

            print('{0:0.0f}% complete\n'.format(perc * 100))
            print('Approximately {0} remaining'
                  .format(format_time_str(remain_time)))
            print(sols)

    # Close out solution timer

    sols.stopwatch(end=True)

    # Save solutions to file if any are found

    if sols.num_sol:
        save_obj(sols, filename=time_str)

    # Delete folder if no ray traceable solutions were found

    if not sols.num_sol_rt:
        folder_exist(path, erase=True)

    return sols


def check_sph_var_sol(efx, oal, bfl, f1, f2, f3_x, f3_y, f4, min_air=0.):

    """


    Identifies whether there is a valid solution for a five group spherical
    variator anamorphic zoom design. Only valid solution types (via Dodoc) are:

        A:     P - N - (Py - Px) - P
        B:     N - N - (Py - Px) - P
        C:     N - P - (Py - Px) - P
        D:     P - P - (Nx - Ny) - P
        E:     N - P - (Nx - Ny) - P

    Only solution types A and E are found for reasonable boundary conditions.

    B needs BFl > 300-ish
    C needs BFL > 200-ish, f3_x > 400-ish, f1 > -20
    D needs BFL < 20-ish, f1 > 400-ish

    I also found a P - N - (Nx - Ny) - P solution that Dodoc didn't include that
    doesn't crash or have internal images...albeit the design was absurdly long

    References:

    [1]  A. Dodoc, “Anamorphic prime and zoom lenses,” in Zoom Lenses VI, San
         Diego, United States, Sep. 2019, p. 2, doi: 10.1117/12.2527911.


    efx:            effective focal lengths in X of zoom system [mm]

    oal:            overall length of zoom system [mm]

    bfl:            back focal length of zoom system [mm]

    f1:             group 1 focal length [mm]

    f2:             group 2 focal length [mm]

    f3_x:           group 3 focal length in X [mm]

    f3_y:           group 3 focal length in Y [mm]

    f4:             group 4 focal length [mm]

    min_air:        minimum air space between groups; default is 0

    """

    # Initialize data structure

    sols = []

    # Calculate group 2 magnification and image conjugates

    m2 = -efx * f3_x / f1 / f4
    s2 = f2 * (1 / m2 - 1)
    s2p = f2 * (1 - m2)

    # Calculate group airspaces

    t1 = f1 - s2
    t2_x = s2p + f3_x
    t2_y = s2p + f3_y

    # Calculate group zoom motion and form focal length and type arrays based
    # on order of cylindrical compensators

    group_z = np.ones((5, efx.size)) * bfl
    group_z[0, :] += oal
    group_z[1, :] = group_z[0, :] - t1

    # YX compensator order (positive f3's)

    if t2_x[0] > t2_y[0]:

        group_z[2, :] = group_z[1, :] - t2_y
        group_z[3, :] = group_z[1, :] - t2_x

        group_efl = np.array([f1, f2, f3_y, f3_x, f4])
        group_type = np.array(['XY', 'XY', 'Y', 'X', 'XY'])

    # XY compensator order (negative f3's)

    else:

        group_z[2, :] = group_z[1, :] - t2_x
        group_z[3, :] = group_z[1, :] - t2_y

        group_efl = np.array([f1, f2, f3_x, f3_y, f4])
        group_type = np.array(['XY', 'XY', 'X', 'Y', 'XY'])

    # Check for a valid solution (no group crashes or internal images). Note
    # that for cylindrical variators there is no easy way to filter internal
    # images.

    if np.diff(group_z[::-1], axis=0).min() > min_air and not 0 < f1 < t1[0]:

        # Add solution to array

        sols.append([group_z, group_efl, group_type])

    return sols
