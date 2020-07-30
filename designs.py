"""
Author:         David Henry Lippman
File:           designs.py
Date created:   07/22/20
Date modified:  07/22/20

"""

from codev import create_ana_zoom, ray_trace, opti_ana_zoom, avg_spot_size, \
    avg_group_efl, avg_clear_aper, tho, save_seq
from utilities import format_time_str
import numpy as np
from time import time
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import plotly.graph_objects as go


class SystemConfig:

    """


    Anamorphic zoom system setup.


    """

    def __init__(self, ana_rat=2., efl_sys_rng=np.array([28., 76.]),
                 bfl_rng=np.array([35., 65.]), ttl_rng=np.array([240., 365.]),
                 efl_group_rng=np.array([20., 500.]), num_zoom=25, fno=4.,
                 img_dim=np.array([22.31, 18.67]),
                 wl=np.array([656.2725, 587.5618, 486.1327]), wl_ref=1):

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

        repr_str = '\n' + 10 * '~' + ' SYSTEM CONFIGURATION ' + 10 * '~' + \
                   '\n\n'

        repr_str += 'Anamorphic ratio:\t\t{0:0.0f}\n'.format(self.ana_rat)
        repr_str += 'EFL, X:\t\t\t\t\t[{0:0.2f}, {1:0.2f}] mm\n' \
            .format(self.efl_sys_rng[0], self.efl_sys_rng[1])
        repr_str += 'EFL, Y:\t\t\t\t\t[{0:0.2f}, {1:0.2f}] mm\n\n' \
            .format(self.ana_rat * self.efl_sys_rng[0],
                    self.ana_rat * self.efl_sys_rng[1])

        repr_str += 'BFL:\t\t\t\t\t[{0:0.2f}, {1:0.2f}] mm\n' \
            .format(self.bfl_rng[0], self.bfl_rng[1])
        repr_str += 'TTL:\t\t\t\t\t[{0:0.2f}, {1:0.2f}] mm\n' \
            .format(self.ttl_rng[0], self.ttl_rng[1])
        repr_str += 'EFL, group:\t\t\t\t[{0:0.2f}, {1:0.2f}] mm\n\n' \
            .format(self.efl_group_rng[0], self.efl_group_rng[1])

        repr_str += 'F/number:\t\t\t\t{0:0.2f}\n'.format(self.fno)
        repr_str += 'Wavelength(s):\t\t\t'
        for w in self.wl:
            repr_str += '{0:0.2f}\t'.format(w)
        repr_str += ' nm\n'
        repr_str += 'HFOV, X x Y:\t\t\t[{0:2.2f} x {1:2.2f},\n' \
                    '\t\t\t\t\t\t {2:2.2f} x {3:2.2f}] deg\n' \
            .format(np.degrees(self.xan[0, -1]),
                    np.degrees(self.yan[0, -1]),
                    np.degrees(self.xan[-1, -1]),
                    np.degrees(self.yan[-1, -1]))

        repr_str += '\n' + (2 * 10 + 22) * '~' + '\n\n'

        return repr_str

    def set_fov(self, scale_fact=1.):

        """


        Sets the system angular field-of-view.


        scale_fact:     scale factor for image space field of view;
                        default is 1

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

        # Loop over first and final zoom positions

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

    def __init__(self, config, sol, num_zoom=5, sol_num=None):

        """


        Initialize anamorphic zoom design object


        config:         system configuration object

        sol:            solution list of lists generated in the Monte Carlo
                        search

        num_zoom:       number of zoom positions to consider for the design

        sol_num:        solution number used as an identifier

        """

        # Initialize variables

        self.config = config
        self.group_z_og = sol[0]
        self.group_efl = sol[1]
        self.group_type = sol[2]
        self.num_zoom = num_zoom
        self.sol_num = sol_num

        if self.num_zoom > self.config.num_zoom:
            self.num_zoom = self.config.num_zoom
            print('WARNING: too few zoom positions for design')

        # Evaluation parameters

        self.ray_trace = 0
        self.avg_spot_size = 0
        self.avg_group_efl = 0
        self.avg_clear_aper = 0
        self.tho = np.empty((5, self.num_zoom))

        # Calculate even spacing indices for number of design zoom positions

        self.ind = np.linspace(0, config.num_zoom - 1, 5).astype(int)

        # Calculate zoom position focal lengths

        self.efx = config.efx[self.ind]
        self.efy = config.efy[self.ind]

        # Calculate zoom positions motion

        self.group_z = self.group_z_og[:, self.ind]

        # Calculate basic values

        self.num_group = self.group_efl.size
        self.ttl = self.group_z[0, 0]
        self.bfl = self.group_z[-1, 0]
        self.oal = self.ttl - self.bfl

        # Classify solution type

        self.classify()

    def __repr__(self):

        repr_str = '\n' + 10 * '~' + ' ANAMORPHIC ZOOM DESIGN ' + 10 * '~' + \
                   '\n\n'

        if self.sol_num is not None:
            repr_str += 'Solution number:\t\t{0:0.0f}\n\n'.format(self.sol_num)

        repr_str += 'EFX:\t\t\t\t\t[{0:0.1f}, {1:0.1f}] mm\n' \
                    .format(self.efx[0], self.efx[-1])
        repr_str += 'EFY:\t\t\t\t\t[{0:0.1f}, {1:0.1f}] mm\n\n' \
                    .format(self.efy[0], self.efy[-1])

        repr_str += 'Variator type:\t\t\t{0}\n'\
                    .format(self.vari_str.capitalize())
        if self.power_type_x == self.power_type_y:
            repr_str += 'Solution type:\t\t\t{0}\n'.format(self.power_type_x)
        else:
            repr_str += 'Solution type:\t\t\t{0}\n'.format(self.power_type)
        repr_str += 'Cyl. orient.:\t\t\t{0}\n\n'.format(self.orient_type)

        repr_str += 'Groups:\n\n'

        template = '{0:>4.0f}{1:>8s}{2:>8s}{3:>14.4f} mm\n'
        for g in range(self.num_group):

            surf_type = 'CYL'
            if self.group_type[g] == 'XY':
                surf_type = 'SPH'

            surf_orient = '-'
            if surf_type == 'CYL':
                surf_orient = self.group_type[g]

            repr_str += template.format(g + 1, surf_type, surf_orient,
                                        self.group_efl[g])

        repr_str += '\n'

        if self.avg_spot_size:
            repr_str += 'Avg. spot size:\t\t\t{0:0.2f} mm\n'\
                        .format(self.avg_spot_size)
        if self.avg_group_efl:
            repr_str += 'Avg. group EFL:\t\t\t{0:0.2f} mm\n'\
                        .format(self.avg_group_efl)
        if self.avg_clear_aper:
            repr_str += 'Avg. clear aper.:\t\t{0:0.2f} mm\n'\
                        .format(self.avg_clear_aper)

        repr_str += '\n' + (2 * 10 + 24) * '~' + '\n\n'

        return repr_str

    def classify(self):

        """


        Classifies the solution type


        """

        # Classify variator type

        if self.num_group == 5:
            self.vari_str = 'spherical'
        elif self.num_group == 6:
            self.vari_str = 'cylindrical'

        # Classify solution and orientation types

        self.power_type = ''
        self.power_type_x = ''
        self.power_type_y = ''
        self.orient_type = ''

        # Loop over groups

        for g in range(self.num_group):

            # Determine group type

            if self.group_type[g] == 'XY':

                # Classify power

                if self.group_efl[g] > 0:
                    self.power_type += 'P'
                    self.power_type_x += 'P'
                    self.power_type_y += 'P'
                else:
                    self.power_type += 'N'
                    self.power_type_x += 'N'
                    self.power_type_y += 'N'

            elif self.group_type[g] == 'X':

                # Classify power

                if self.group_efl[g] > 0:
                    self.power_type += 'P'
                    self.power_type_x += 'P'
                else:
                    self.power_type += 'N'
                    self.power_type_x += 'N'

                # Classify cylinder orientation

                self.orient_type += self.group_type[g]

            elif self.group_type[g] == 'Y':

                # Classify power

                if self.group_efl[g] > 0:
                    self.power_type += 'P'
                    self.power_type_y += 'P'
                else:
                    self.power_type += 'N'
                    self.power_type_y += 'N'

                # Classify cylinder orientation

                self.orient_type += self.group_type[g]

        # Create full solution type string

        self.sol_type = self.power_type + '_' + self.orient_type

    def make_codev(self, save_seq=None):

        """


        Creates anamorphic zoom design in CODE V


        save_seq:       filename to save the design as a sequence file

        """

        create_ana_zoom(self, save_seq)

    def check_ray_trace(self):

        """


        Checks for ray trace errors in CODE V model


        """

        # Checks rays and stores result

        self.ray_trace = ray_trace()

        return self.ray_trace

    def optimize(self):

        """


        Optimizes an anamorphic zoom design in CODE V


        """

        opti_ana_zoom()

    def analyze(self):

        """


        Analyze CODE V designs for performance, aberrations, and packaging


        """

        # Calculate average spot size

        self.avg_spot_size = avg_spot_size()

        # Calculate average group efl

        self.avg_group_efl = avg_group_efl(self.num_group)

        # Calculate average group efl

        self.avg_clear_aper = avg_clear_aper()

        # Calculate third order aberration values

        self.tho = tho()

    def save_seq(self, path=None, out=False):

        """


        Save anamorphic zoom design as CODE V sequence file


        filename:       filename to save

        """

        # Check if path was provided

        if path is None:
            path = "C:\\CVUSER\\"

        # Create filename if not provided

        filename = '{0}_sol{1}.seq'.format(self.sol_type, self.sol_num)

        # Create full path

        path += filename

        # Save sequence file

        save_seq(path)

        # Print message, if desired

        if out:
            print('Solution saved:\t{0}'.format(path))

    def plot_zoom(self):

        """


        Plot zoom motion of design


        """

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

            ax.plot(self.group_z_og[i, :], self.config.efx, color=clrs[i],
                    linestyle=lstyles[i], label=leg_str)

            ax_par.plot(self.group_z_og[i, :], self.config.efy, color='none')

        # Set plot settings

        ax.set(title=self.sol_type,
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


class Solutions:

    """


    Anamorphic zoom solutions object


    """

    def __init__(self, config, num_trial, var_type, same_xy):

        """

        Initialize anamorphic zoom solutions object


        """

        # Initialize variables

        self.config = config
        self.num_trial = num_trial
        self.var_type = var_type[0:3].upper()
        self.same_xy = same_xy

        self.sols = []
        self.sols_rt = []

        self.num_sol = 0
        self.num_sol_rt = 0

        # Create lists of all possible power and orientation types based on
        # variator type

        if self.var_type == "CYL":
            power_list = [c1 + c2 + c3 + c4 + c5 + c6 for c1 in ["P", "N"]
                                                      for c2 in ["P", "N"]
                                                      for c3 in ["P", "N"]
                                                      for c4 in ["P", "N"]
                                                      for c5 in ["P", "N"]
                                                      for c6 in ["P", "N"]]
            orient_list = ["XXYY", "XYXY", "XYYX", "YXXY", "YXYX", "YYXX"]

        elif self.var_type == "SPH":
            power_list = [c1 + c2 + c3 + c4 + c5 for c1 in ["P", "N"]
                                                 for c2 in ["P", "N"]
                                                 for c3 in ["P", "N"]
                                                 for c4 in ["P", "N"]
                                                 for c5 in ["P", "N"]]
            orient_list = ["XY", "YX"]

        else:
            sys.exit('Invalid variator type')

        # Create empty dictionary of all possible solution types

        self.power_type = {}
        self.power_type_rt = {}
        self.orient_type = {}
        self.orient_type_rt = {}
        self.sol_type = {}
        self.sol_type_rt = {}

        # Loop over power types

        for power_key in power_list:

            self.power_type[power_key] = 0
            self.power_type_rt[power_key] = 0

            # Loop over orientation types

            for orient_key in orient_list:

                self.orient_type[orient_key] = 0
                self.orient_type_rt[orient_key] = 0

                self.sol_type[power_key + "_" + orient_key] = 0
                self.sol_type_rt[power_key + "_" + orient_key] = 0

        # Start stopwatch

        self.start_time = time()
        self.stop_time = self.start_time
        self.comp_time = self.stop_time - self.start_time
        self.sw_status = True

    def __repr__(self):

        # Update stop watch, if still active

        self.stopwatch()

        # Header

        repr_str = '\n' + 10 * '~' + ' ANAMORPHIC ZOOM SOLUTIONS ' + \
                   10 * '~' + '\n\n'

        # First order solutions, if any found

        if self.num_sol:
            repr_str += 'First order solutions:\t{0:0.0f}\n\n'\
                        .format(self.num_sol)
            template = '{0:>12s}{1:>8.0f}\n'
            for key in self.sort_sols(self.sol_type):
                if self.sol_type[key]:
                    power_type, orient_type = self.split_sol_type(key)
                    power_type = self.abbrev_power_type(power_type, orient_type)
                    sol_str = power_type + ' ' + orient_type
                    repr_str += template.format(sol_str, self.sol_type[key])

            # Ray traceable solutions order solutions, if any found

            if self.num_sol_rt:
                repr_str += '\nRay traceable solutions:\t{0:0.0f}\n\n'\
                            .format(self.num_sol_rt)
                for key in self.sort_sols(self.sol_type_rt):
                    if self.sol_type_rt[key]:
                        power_type, orient_type = self.split_sol_type(key)
                        power_type = self.abbrev_power_type(power_type,
                                                            orient_type)
                        sol_str = power_type + ' ' + orient_type
                        repr_str += template.format(sol_str,
                                                    self.sol_type_rt[key])

        # No solutions found

        else:
            repr_str += 'No solutions found\n'


        # Print computation time

        repr_str += '\nComputation time:\t\t{0}\n'\
                    .format(format_time_str(self.comp_time))

        repr_str += '\n' + (2 * 10 + 27) * '~' + '\n\n'

        return repr_str

    def add_sol(self, sol):

        """


        Adds solution to solution set


        sol:        anamorphic zoom solution to add

        """

        # Add general solution

        self.num_sol += 1
        self.sols.append(sol)
        self.power_type[sol.power_type] += 1
        self.orient_type[sol.orient_type] += 1
        self.sol_type[sol.sol_type] += 1

        # Add ray traceable solution, if applicable

        if sol.ray_trace:

            self.num_sol_rt += 1
            self.sols_rt.append(sol)
            self.power_type_rt[sol.power_type] += 1
            self.orient_type_rt[sol.orient_type] += 1
            self.sol_type_rt[sol.sol_type] += 1

    def get_sol_num(self, sol_num):

        """


        Gets a solutions based on solution index. This solution number is
        zero-indexed.


        sol_num:        solution number zero-indexed


        """

        return self.sols[sol_num]

    def sort_sols(self, sols):

        """


        Sorts a provided solution dictionary by value, descending


        sols:       solution dictionary to sort

        """

        return {k: v for k, v in sorted(sols.items(),
                                        key=lambda item: item[1],
                                        reverse=True)}

    def split_sol_type(self, sol_type):

        """


        Splits a solution type string into power and orientation strings.


        sol_type:       solution type string to split

        """

        # Cylindrical variator

        if self.var_type == "CYL":
            return sol_type[0: 6], sol_type[-4:]

        # Spherical variator

        return sol_type[0: 5], sol_type[-2:]

    def abbrev_power_type(self, power_type, orient_type):

        """


        Abbreviates the power type string if the solution type is the same in X
        and Y


        power_type:         power type string

        orient_type:        orientation type string

        """

        # If not the same solution type in X and Y, no abbreviation can be
        # performed

        if not self.same_xy:
            return power_type

        # If the same solution type in X and Y, abbreviation can be performed

        new_power_type = power_type[0]
        if self.var_type == "SPH":
            new_power_type += power_type[1]
        ind = len(new_power_type)
        for i, c in enumerate(orient_type):
            if c == "X":  # combine X and Y power types based off of X
                new_power_type += power_type[ind + i]
        new_power_type += power_type[-1]

        return new_power_type

    def stopwatch(self, end=False):

        """


        Updates the stopwatch on the execution duration


        """

        # Update computation time

        if self.sw_status:
            self.stop_time = time()
            self.comp_time = self.stop_time - self.start_time

        # Stop stopwatch

        if end:
            self.sw_status = False

        return self.comp_time

    def demographics(self):

        # Initialize data structures

        sol_type = []
        sol_type_rt = []
        ray_trace = {}
        avg_spot_size = {}
        avg_group_efl = {}
        avg_clear_aper = {}
        SA = {}
        TCO = {}
        TAS = {}
        SAS = {}
        PTB = {}

        # Determine ray traceability likelihood for found first order solutions

        for sol in self.sols:

            # Check to see if this solution type has been initialized yet

            if sol.sol_type not in sol_type:

                sol_type.append(sol.sol_type)
                ray_trace[sol.sol_type] = np.zeros(2, dtype=int)

            # Check whether solution ray traces

            if sol.ray_trace:
                ray_trace[sol.sol_type][0] += 1
            else:
                ray_trace[sol.sol_type][1] += 1

        # Calculate ray traceability liklihood

        for s_type in sol_type:
            ray_trace[s_type] = 100 * ray_trace[s_type][0] / \
                                      ray_trace[s_type].sum()

        # Loop over ray traceable solutions

        for sol in self.sols_rt:

            # Check to see if this solution type has been initialized yet

            if sol.sol_type not in sol_type_rt:

                # Add evaluation parameters to dictionary by initializing list

                sol_type_rt.append(sol.sol_type)
                avg_spot_size[sol.sol_type] = [sol.avg_spot_size]
                avg_group_efl[sol.sol_type] = [sol.avg_group_efl]
                avg_clear_aper[sol.sol_type] = [sol.avg_clear_aper]
                SA[sol.sol_type] = [sol.tho[0, :].mean()]
                TCO[sol.sol_type] = [sol.tho[1, :].mean()]
                TAS[sol.sol_type] = [sol.tho[2, :].mean()]
                SAS[sol.sol_type] = [sol.tho[3, :].mean()]
                PTB[sol.sol_type] = [sol.tho[4, :].mean()]

            else:

                # Add evaluation parameters to dictionary by appending to list

                avg_spot_size[sol.sol_type].append(sol.avg_spot_size)
                avg_group_efl[sol.sol_type].append(sol.avg_group_efl)
                avg_clear_aper[sol.sol_type].append(sol.avg_clear_aper)
                SA[sol.sol_type].append(sol.tho[0, :].mean())
                TCO[sol.sol_type].append(sol.tho[1, :].mean())
                TAS[sol.sol_type].append(sol.tho[2, :].mean())
                SAS[sol.sol_type].append(sol.tho[3, :].mean())
                PTB[sol.sol_type].append(sol.tho[4, :].mean())

        # Plot ray traceability by solution type

        fig, ax = plt.subplots()

        ax.bar(ray_trace.keys(), ray_trace.values())
        ax.set(xlabel='Solution type', ylabel='Ray traceability likelihood [%]')

        # Plot average spot size and average group EFL by ray traceable solution
        # type

        fig, ax = plt.subplots()

        for s_type in sol_type_rt:
            ax.scatter(avg_spot_size[s_type], avg_group_efl[s_type], s=5,
                       label=s_type)
            ax.set(xlabel='Average spot size [mm]',
                   ylabel='Average group EFL [mm]')
            ax.set_xlim(0, 5)
            ax.legend()

        """
        
        Others:
        
        - cylinder orientation
        - aberration breakdown
        
        """

        return

    def plot_sankey(self):

        # Initialize variables

        plot_fail_factor = 5  # how many times more unfound solutions to plot
                              # not to scale

        # Determine number of trials string

        if self.num_trial >= 1e9:
            num_trial_str = '{0:0.0f} billion'.format(self.num_trial / 1e9)
        elif self.num_trial >= 1e6:
            num_trial_str = '{0:0.0f} million'.format(self.num_trial / 1e6)
        elif self.num_trial >= 1e3:
            num_trial_str = '{0:0.0f} thousand'.format(self.num_trial / 1e3)
        else:
            num_trial_str = '{0:0.0f} '.format(self.num_trial)

        # Abbreviate power type, if possible, for plotting clarity

        # First order solutions

        power_type = {}

        for key, value in zip(self.sol_type.keys(), self.sol_type.values()):

            # Form abbreviated power label

            power_str, orient_str = self.split_sol_type(key)
            power_str_abbrev = self.abbrev_power_type(power_str, orient_str)

            # Add/update dictionary power label
            try:
                power_type[power_str_abbrev] += value
            except KeyError:
                power_type[power_str_abbrev] = value

        # Ray traceable solutions

        power_type_rt = {}

        for key, value in zip(self.sol_type_rt.keys(),
                              self.sol_type_rt.values()):

            # Form abbreviated power label

            power_str, orient_str = self.split_sol_type(key)
            power_str_abbrev = self.abbrev_power_type(power_str, orient_str)

            # Add/update dictionary power label
            try:
                power_type_rt[power_str_abbrev] += value
            except KeyError:
                power_type_rt[power_str_abbrev] = value

        # Sort solution types by number found, for plotting clarity

        power_type_sort = self.sort_sols(power_type)
        power_type_rt_sort = self.sort_sols(power_type_rt)

        # Form data arrays for Sankey plot

        label = ["{0} Monte Carlo trials".format(num_trial_str),
                 "Invalid solutions", "Invalid solutions"]
        color = ["blue", "red", "red"]

        source = [0, 1]
        target = [1, 2]
        value = [0, 0]

        # Loop over solutions

        num_sol_fo_count = 0
        node_count = 2  # zero-indexed

        for key, val in zip(power_type_sort.keys(), power_type_sort.values()):

            # Initialize variables

            num_sol_fo = val
            num_sol_rt = power_type_rt_sort[key]
            num_sol_fo_rt_fail = num_sol_fo - num_sol_rt

            # Add first order solution node

            if num_sol_fo:
                label.append("{0} solutions ({1:0.0f})".format(key, num_sol_fo))
                color.append("green")

                node_count += 1
                num_sol_fo_count += num_sol_fo

                source.append(0)
                target.append(node_count)
                value.append(num_sol_fo)

            # Add first order solutions that don't ray trace

            if num_sol_fo_rt_fail:
                source.append(node_count)
                target.append(2)
                value.append(num_sol_fo_rt_fail)

            # Add ray traceable solution node

            if num_sol_rt:
                label.append("{0} solutions ({1:0.0f})".format(key, num_sol_rt))
                color.append("green")

                source.append(node_count)
                node_count += 1
                target.append(node_count)
                value.append(num_sol_rt)

        # Add first value for invalid first order solutions based on fail factor
        # for not to scale plotting purposes

        value[0] = plot_fail_factor * num_sol_fo_count
        value[1] = plot_fail_factor * num_sol_fo_count

        # Form node and link dictionaries

        node = dict(pad=25,
                    thickness=30,
                    line=dict(color=None),
                    label=label,
                    color=color)

        link = dict(source=source,
                    target=target,
                    value=value)

        # Create Sankey data

        data = [go.Sankey(node=node, link=link)]

        # Initialize Sankey figure

        fig = go.Figure(data=data)

        # Set plot settings

        title = "Monte Carlo First-Order Search"
        title += " - {0} Variator".format(self.sols[0].vari_str.capitalize())

        font_size = 10
        fig.update_layout(title_text=title,
                          font=dict(family="Arial", size=font_size))
        fig.show()

        

        return


