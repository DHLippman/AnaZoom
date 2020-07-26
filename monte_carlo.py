"""
Author:         David Henry Lippman
File:           monte_carlo.py
Date created:   07/22/20
Date modified:  07/22/20

"""

from designs import AnamorphicZoom, Solutions
from utilities import rand_rng, get_time_str, folder_exist, save_obj
import numpy as np
from time import time


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

    # Create solutions object

    sols = Solutions(config, num_trial)

    # Create folder to store CODE V solutions in

    time_str = get_time_str()
    path = 'C:\\CVUSER\\Anamorphic Zoom Solutions\\Cylindrical Variator\\' \
           + time_str + '\\'
    folder_exist(path)

    # Set up progress bar and timer

    progress = (num_trial * np.linspace(0.1, 1, 10)).astype(int)
    start_time = time()

    # Loop over Monte Carlo trials

    print('\nPerforming a Monte Carlo search with {0:0.1e} trials\n'
          .format(num_trial))

    for i in range(num_trial):

        # Randomly pick a design form

        design_form = design_forms[np.random.randint(0, design_forms.shape[0])]

        # Randomly pick values for system TTL and BFL

        ttl = rand_rng(config.ttl_rng[0], config.ttl_rng[1], sign=1)
        bfl = rand_rng(config.bfl_rng[0], config.bfl_rng[1], sign=1)
        oal = ttl - bfl

        # Look for solution(s) in X

        # Randomly pick group focal lengths for the stationary groups f1 and f4
        # and one of the two moving groups f2

        f1 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[0])
        f2_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[1])
        f3_x = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                        sign=design_form[2])
        f4 = rand_rng(config.efl_group_rng[0], config.efl_group_rng[1],
                      sign=design_form[3])

        # Attempt to randomly pick a value for the group focal length of the
        # second moving group f3 while ensuring that a zoom solution will be
        # found

        #f3_x = rand_f3(config.efx, oal, bfl, f1, f2_x, f4, config.efl_group_rng,
        #               design_form)

        # Check for valid four group zoom solution(s), if possible

        sols_x = check_four_group(config.efx, oal, bfl, f1, f2_x, f3_x, f4,
                                  min_air=min_air)

        # Loop over solutions in X, if any were found

        for sol_x in sols_x:

            # Look for a compatible Y solution

            sol_y, \
            f2_y, \
            f3_y, = find_y_sol(config, design_form, sol_x, oal, bfl, f1, f4,
                               min_air, samp=10)

            # Combine X solution and Y solution (if found)

            if sol_y is not None:

                """
                
                Successful first order solution
                
                """

                # Combine X and Y solutions for group focal lengths, zoom
                # motions, and group types

                group_efl, \
                group_z, \
                group_type = combine_sols(sol_x, sol_y, f1, f2_x, f2_y, f3_x,
                                          f3_y, f4)

                # Create anamorphic zoom design

                ana_zoom = AnamorphicZoom(config, group_efl, group_type,
                                          group_z)

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

                        ana_zoom.save_seq(sol_num=sols.num_sol, path=path)

                # Add design to solutions

                sols.add_sol(ana_zoom)

        # Check status

        if i + 1 in progress:

            perc = (np.asscalar(np.argwhere(i + 1 == progress)) + 1) / \
                   progress.size
            cur_time = time()
            pass_time = cur_time - start_time
            full_time = pass_time / perc
            remain_time = full_time - pass_time
            print('{0:0.0f}% complete\n'.format(perc * 100))
            if remain_time > 3600:
                print('Approximately {0:0.2f} hours remaining'
                      .format(remain_time / 3600))
            elif remain_time > 60:
                print('Approximately {0:0.2f} minutes remaining'
                      .format(remain_time / 60))
            elif remain_time > 0:
                print('Approximately {0:0.2f} seconds remaining'
                      .format(remain_time))
            print(sols)

    # Save solutions to file if any are found

    if sols.num_sol:
        save_obj(sols, filename=time_str)

    # Delete folder if no ray traceable solutions were found

    if not sols.num_sol_rt:
        folder_exist(path, erase=True)

    return sols


def rand_f3(efl, oal, bfl, f1, f2, f4, efl_group_rng, design_form):

    # Initialize variables

    efl_group_rng_signed = efl_group_rng * design_form[2]  # for group 3
    f3_min = efl_group_rng_signed.min()
    f3_max = efl_group_rng_signed.max()

    # Calculate quadratic equation coefficients

    s4 = -f4 * bfl / (bfl - f4)
    L = oal - f1 + s4
    M = -efl * (oal - L - f1) / (bfl * f1)

    # Calculate maximum allowable value of f3 for which a solution will be found

    term_1 = L ** 2 / 4 - L * f2
    term_2 = L + f2 * (M - 1) ** 2 / M

    if term_2.min() > 0:
        if np.min(term_1 / term_2) < f3_max:
            f3_max = np.min(term_1 / term_2)

    elif term_2.max() < 0:
        if np.max(term_1 / term_2) > f3_min:
            f3_min = np.max(term_1 / term_2)

    else:
        return None

    # Randomly choose an allowable value of f3, if possible

    if f3_min > f3_max:
        return None

    return rand_rng(f3_min, f3_max, sign=design_form[2])



def check_six_group(efx, efy, oal, bfl, f1, f2_x, f2_y, f3_x, f3_y, f4,
                    min_air=0., sol_check=None):

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
    M_x = -efx * (oal - L - f1) / (bfl * f1)
    c_x = L * (f2_x + f3_x) + ((M_x - 1) ** 2) * f2_x * f3_x / M_x
    M_y = -efy * (oal - L - f1) / (bfl * f1)
    c_y = L * (f2_y + f3_y) + ((M_y - 1) ** 2) * f2_y * f3_y / M_y

    # Check if roots are real for all zoom positions (all values of c) for both
    # the X and Y slices

    if L ** 2 - 4 * c_x.min() >= 0 and L ** 2 - 4 * c_x.max() >= 0 and \
       L ** 2 - 4 * c_y.min() >= 0 and L ** 2 - 4 * c_y.max() >= 0:

        # Check positive root in X

        t2_x = (L + np.sqrt(L ** 2 - 4 * c_x)) / 2

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


def calc_zoom_motion(efx, efy, oal, bfl, M_x, M_y, t2_x, t2_y, f2_x, f2_y, L,
                     s4):

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

    # Calculate air spaces in X slice

    s2_x = ((M_x - 1) * t2_x + L) / ((M_x - 1) - M_x * t2_x / f2_x)
    t1_x = oal - L + s4 - s2_x
    t3_x = oal - t1_x - t2_x

    # Calculate air spaces in Y slice

    s2_y = ((M_y - 1) * t2_y + L) / ((M_y - 1) - M_y * t2_y / f2_y)
    t1_y = oal - L + s4 - s2_y
    t3_y = oal - t1_y - t2_y

    # Calculate group zoom motions

    z = np.ones((6, efx.size)) * bfl
    z[4, :] += t3_x
    z[3, :] += t3_y
    z[2, :] = t2_x + z[4, :]
    z[1, :] = t2_y + z[3, :]
    z[0, :] += oal

    # Sort moving group ordering based on arbitrary (first) zoom position

    z_moving = z[1: -1, :]
    z[1: -1, :] = z_moving[z_moving[:, 0].argsort()[::-1]]

    return z


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


def calc_zoom_motion_old(M, L, t2, oal, bfl, f2, s4, efl):

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


def find_y_sol(config, design_form, sol_x, oal, bfl, f1, f4, min_air,
               samp=10):

    """


    Finds Y solutions that are compatible with a given X solution.

    Routine has the potential to find multiple compatible Y solutions for the
    given X solution, but in order to not introduce bias (multiple full
    solutions with the same X group focal lengths but different Y ones) only the
    first found solution is accepted. For the same reason, the ordering of the
    potential f2, f3 solutions are shuffled


    Increases the likelihood of finding a compatible Y solution by approximately
    10 x

    """

    # Find group EFL bounds to randomly check between

    bounds = np.linspace(config.efl_group_rng.min(),
                         config.efl_group_rng.max(), samp + 1)

    # Determine random f2, f3 values to check

    f2_check = np.empty(samp)
    f3_check = np.empty(samp)

    for i in range(samp):
        f2_check[i] = rand_rng(bounds[i], bounds[i + 1], sign=design_form[1])
        f3_check[i] = rand_rng(bounds[i], bounds[i + 1], sign=design_form[2])

    # Flatten a mesh into a 2 column array with f2 and f3 group focal lengths

    f2_mesh, f3_mesh = np.meshgrid(f2_check, f3_check)

    f_search = np.hstack((np.reshape(f2_mesh, (samp ** 2, 1)),
                          np.reshape(f3_mesh, (samp ** 2, 1))))

    # Randomly shuffle group focal lengths to search over

    np.random.shuffle(f_search)

    # Loop over randomly-ordered f2, f3 combinations looking for a compatible
    # solution

    for f2_y, f3_y in f_search:

        # Check for valid four group zoom solution(s) taking into
        # account solution found in X

        sols_y = check_four_group(config.efy, oal, bfl, f1, f2_y, f3_y,
                                  f4, min_air=min_air, sol_check=sol_x)

        # Valid solution(s) found

        for sol_y in sols_y:

            # Return the first of potential two solutions (for +/- roots)

            return sol_y, f2_y, f3_y

    # If no solution is found, return None instead

    return None, None, None


def combine_sols(sol_x, sol_y, f1, f2_x, f2_y, f3_x, f3_y, f4):

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

