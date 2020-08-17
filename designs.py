"""
Author:         David Henry Lippman
File:           designs.py
Date created:   07/22/20
Date modified:  08/17/20

"""

from codev import create_ana_zoom, ray_trace, opti_ana_zoom, avg_spot_size, \
    avg_group_efl, avg_clear_aper, tho, save_seq
from utilities import format_time_str, sort_dict
import numpy as np
from time import time
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
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

        # Calculate OAL range given TTL and BFL range

        self.oal_rng = self.ttl_rng - np.flip(self.bfl_rng, axis=0)

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

        self.ray_trace = False
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

        repr_str += 'TTL:\t\t\t\t\t{0:0.2f} mm\n'.format(self.ttl)
        repr_str += 'BFL:\t\t\t\t\t{0:0.2f} mm\n\n'.format(self.bfl)

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
        repr_str += 'Ray traceable:\t\t\t{0}\n'.format(self.ray_trace)

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

        # Calculates average spot size across all fields and zooms

        self.avg_spot_size = avg_spot_size()

        # Calculates average group efl across all groups

        self.avg_group_efl = avg_group_efl(self.num_group)

        # Calculates average clear aperture across all surfaces and zooms

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
        lstyles[self.group_type == 'X'] = ':'
        lstyles[self.group_type == 'Y'] = '--'

        # Initialize figure

        font = {'size': 14}
        plt.rc('font', **font)
        fig, ax = plt.subplots(figsize=(12, 6))
        plt.subplots_adjust(left=0.25, right=0.65)

        # Create second y-axis

        ax_par = ax.twinx()
        ax_par.spines["left"].set_position(("axes", -0.15))
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

            if self.vari_str == "CYL":
                group_ind = np.ceil(i / 2).astype(int) + 1
            else:
                if i < 3:
                    group_ind = i + 1
                else:
                    group_ind = i

            type_str = self.group_type[i]
            if type_str == 'XY':
                type_str = ''

            subscript_str = '{0}'.format(group_ind) + type_str
            leg_str = r'$f_{{{0:2s}}}$ = {1:0.2f} mm'.format(subscript_str,
                                                             self.group_efl[i])

            # Plot zoom motion

            ax.plot(self.group_z_og[i, :], self.config.efx, color=clrs[i],
                    linestyle=lstyles[i], label=leg_str, lw=3)

            ax_par.plot(self.group_z_og[i, :], self.config.efy, color='none',
                        lw=3)

        # Set plot settings

        ax.set(title=self.sol_type,
               xlabel="Group position, z [mm]",
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


    Anamorphic zoom solutions object. Note that with plotly renderer, nodes can
    be manually re-positioned.


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

        # Solutions info  

        repr_str += 'Number of trials:\t\t{0:0.2e}\n'.format(self.num_trial)
        repr_str += 'Variator type:\t\t\t{0}\n'.format(self.var_type)
        repr_str += 'Same X-Y type:\t\t\t{0}\n'.format(self.same_xy)
        repr_str += 'Computation time:\t\t{0}\n\n'\
                    .format(format_time_str(self.comp_time))

        # First order solutions, if any found

        if self.num_sol:
            repr_str += 'First order solutions:\t{0:0.0f}\n\n'\
                        .format(self.num_sol)
            template = '{0:>12s}{1:>8.0f}\n'
            for key in sort_dict(self.sol_type):
                if self.sol_type[key]:
                    power_type, orient_type = self.split_sol_type(key)
                    power_type = self.abbrev_power_type(power_type, orient_type)
                    sol_str = power_type + ' ' + orient_type
                    repr_str += template.format(sol_str, self.sol_type[key])

            # Ray traceable solutions order solutions, if any found

            if self.num_sol_rt:
                repr_str += '\nRay traceable solutions:\t{0:0.0f}\n\n'\
                            .format(self.num_sol_rt)
                for key in sort_dict(self.sol_type_rt):
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

    def type_sankey(self):

        """


        Makes a Sankey diagram of the Monte Carlo search process by solution
        space type.


        """

        # Initialize variables

        plot_fail_factor = 1  # how many times more unfound solutions to plot
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

        # Remove PPPP due to internal image

        if "PPPP" in power_type:
            power_type.pop("PPPP")

        if "PPPP" in power_type_rt:
            power_type_rt.pop("PPPP")

        # Sort solution types by number found, for plotting clarity

        power_type_sort = sort_dict(power_type)
        power_type_rt_sort = sort_dict(power_type_rt)

        # Form data arrays for Sankey plot

        label = ["{0} Monte Carlo trials".format(num_trial_str),
                 "Invalid solutions", "Invalid solutions"]
        color = ["black", "red", "red"]

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
                color.append("blue")

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
                    thickness=50,
                    line=dict(color=None),
                    label=label,
                    color=color)

        link = dict(source=source,
                    target=target,
                    value=value)

        # Create Sankey data

        data = [go.Sankey(node=node, link=link, arrangement='snap')]

        # Initialize Sankey figure

        fig = go.Figure(data=data)

        # Set plot settings

        title = "Monte Carlo First-Order Search"
        title += " - {0} Variator".format(self.sols[0].vari_str.capitalize())
        title = ''

        font_size = 30
        fig.update_layout(title_text=title,
                          font=dict(family="Helvetica", size=font_size))
        fig.show()

        return

    def type_breakdown(self):

        """


        Plots breakdown of solution types based on power and cylinder
        orientation


        """

        # Create orientation-to-index dictionary

        orient_ind = {'XXYY': 0,
                      'XYXY': 1,
                      'XYYX': 2,
                      'YYXX': 3,
                      'YXYX': 4,
                      'YXXY': 5}
        
        # Create orientation hatching patterns

        density = 3
        orient_hatch = [density * '*', density * '.',
                        density * "/", density * "\\",
                        density * "x", density * "+"]


        # Form dictionaries of valid first-order solutions for power and
        # orientation

        power_dict = {}

        for key, val in zip(self.sol_type.keys(), self.sol_type.values()):

            # Valid first-order solution type found

            if val:

                # Form orientation and power keys

                orient_key = key[-4:]
                if self.var_type == "CYL":
                    if self.same_xy:
                        power_key = key[0]
                        for i, c in enumerate(orient_key):
                            if c == 'X':
                                power_key += key[i + 1]
                        power_key += key[5]
                    else:
                        power_key = key[0: 6]

                else:
                    power_key = key[0: 3] + key[4]

                # Update power dictionary

                try:
                    power_dict[power_key][0, orient_ind[orient_key]] += val
                except KeyError:
                    power_dict[power_key] = np.zeros((2, 6))
                    power_dict[power_key][0, orient_ind[orient_key]] += val

        # Add ray traceable solutions to dictionaries

        for key, val in zip(self.sol_type_rt.keys(), self.sol_type_rt.values()):

            # Ray traceable solution type found

            if val:

                # Form orientation and power keys

                orient_key = key[-4:]
                if self.var_type == "CYL":
                    if self.same_xy:
                        power_key = key[0]
                        for i, c in enumerate(orient_key):
                            if c == 'X':
                                power_key += key[i + 1]
                        power_key += key[5]
                    else:
                        power_key = key[0: 6]

                else:
                    power_key = key[0: 3] + key[4]

                # Update power dictionary

                power_dict[power_key][1, orient_ind[orient_key]] += val

        # Remove PPPP due to internal image

        if "PPPP" in power_dict:
            power_dict.pop("PPPP")

        # Determine dictionary sort order based no the number of valid first-
        # order solutions

        sol_types = np.array(list(power_dict.keys()), dtype='<U6')
        num_sol_fo = np.array([val[0, :].sum()
                               for val in list(power_dict.values())], dtype=int)
        sol_types_sort = sol_types[num_sol_fo.argsort()[::-1]]

        # Make solution array to make stacked bar chart

        num_sol_types = len(sol_types)
        num_sols_fo = np.zeros((6, num_sol_types), dtype=int)
        num_sols_rt = np.zeros((6, num_sol_types), dtype=int)

        for i, sol_type in enumerate(sol_types_sort):
            num_sols_fo[:, i] = power_dict[sol_type][0, :].transpose()
            num_sols_rt[:, i] = power_dict[sol_type][1, :].transpose()

        # Initialize figure

        font = {'size': 14}
        plt.rc('font', **font)
        fig, ax1 = plt.subplots(figsize=(12, 6))
        ax2 = ax1.twinx()

        color1 = '#1b75bb'  # blue
        color2 = '#37b34a'  # green

        ind = np.arange(num_sol_types)  # the label locations
        width = 0.4  # the width of the bars

        # Plot bar charts

        # Loop over orientation rows

        for i, (orient_fo, orient_rt) in enumerate(zip(num_sols_fo,
                                                       num_sols_rt)):

            # Valid first-order solutions

            ax1.bar(ind - width / 2, orient_fo, width=width,
                    bottom=np.sum(num_sols_fo[0: i, :], axis=0),
                    facecolor=color1, edgecolor='k', hatch=orient_hatch[i])

            # Ray traceable solutions

            ax2.bar(ind + width / 2, orient_rt, width=width,
                    bottom=np.sum(num_sols_rt[0: i, :], axis=0),
                    facecolor=color2, edgecolor='k', hatch=orient_hatch[i])

        # Plot settings

        ax1.set_ylabel('Valid first-order solutions', color=color1)
        ax2.set_ylabel('Ray traceable solutions', color=color2)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax2.tick_params(axis='y', labelcolor=color2)
        ax1.set_ylim(0, 1.08 * np.sum(num_sols_fo, axis=0).max())
        ax2.set_ylim(0, 1.08 * np.sum(num_sols_rt, axis=0).max())

        ax1.set_xticks(ind)
        ax1.set_xticklabels(sol_types_sort, rotation=45)

        # Make custom legend

        legend_elements = []
        for hatch, key in zip(orient_hatch, orient_ind.keys()):
            legend_elements.append(ptc.Patch(facecolor='w', edgecolor='k',
                                             hatch=hatch, label=key))
        legend_elements = legend_elements[::-1]

        ax1.legend(handles=legend_elements)

        return

    def type_breakdown_sep(self):

        """


        Plots breakdown of solution types based on power and cylinder
        orientation in separate bar charts


        """

        # Form dictionaries of valid first-order solutions for power and orientation

        power_dict = {}
        orient_dict = {}

        for key, val in zip(self.sol_type.keys(), self.sol_type.values()):

            # Valid first-order solution type found

            if val:

                # Form orientation and power keys

                orient_key = key[-4:]
                if self.var_type == "CYL":
                    if self.same_xy:
                        power_key = key[0]
                        for i, c in enumerate(orient_key):
                            if c == 'X':
                                power_key += key[i + 1]
                        power_key += key[5]
                    else:
                        power_key = key[0: 6]

                else:
                    power_key = key[0: 3] + key[4]

                # Update power dictionary

                try:
                    power_dict[power_key] += val
                except KeyError:
                    power_dict[power_key] = val

                # Update orientation dictionary

                try:
                    orient_dict[orient_key] += val
                except KeyError:
                    orient_dict[orient_key] = val

        # Sort dictionaries by value, descending

        power_dict = sort_dict(power_dict)
        orient_dict = sort_dict(orient_dict)

        # Expand dictionaries to have list values for ray tracing results also

        for key in power_dict:
            power_dict[key] = [power_dict[key], 0]

        for key in orient_dict:
            orient_dict[key] = [orient_dict[key], 0]

        # Add ray traceable solutions to dictionaries

        for key, val in zip(self.sol_type_rt.keys(), self.sol_type_rt.values()):

            # Ray traceable solution type found

            if val:

                # Form orientation and power keys

                orient_key = key[-4:]
                if self.var_type == "CYL":
                    if self.same_xy:
                        power_key = key[0]
                        for i, c in enumerate(orient_key):
                            if c == 'X':
                                power_key += key[i + 1]
                        power_key += key[5]
                    else:
                        power_key = key[0: 6]

                else:
                    power_key = key[0: 3] + key[4]

                # Update power dictionary

                power_dict[power_key][1] += val

                # Update orientation dictionary

                orient_dict[orient_key][1] += val

        # Remove PPPP due to internal image

        if "PPPP" in power_dict:
            power_dict.pop("PPPP")

        # Form power lists for bar chart

        power_labels = []
        sols_power_fo = []
        sols_power_rt = []

        for key, val in zip(power_dict.keys(), power_dict.values()):
            power_labels.append(key)
            sols_power_fo.append(val[0])
            sols_power_rt.append(val[1])

        # Form orientation lists for bar chart

        orient_labels = []
        sols_orient_fo = []
        sols_orient_rt = []

        for key, val in zip(orient_dict.keys(), orient_dict.values()):
            orient_labels.append(key)
            sols_orient_fo.append(val[0])
            sols_orient_rt.append(val[1])

        # Initialize figure

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.5))
        plt.subplots_adjust(hspace=0.4)
        ax1b = ax1.twinx()
        ax2b = ax2.twinx()

        color1 = '#1b75bb'
        color2 = '#37b34a'

        # Plot bar charts

        # Power

        x = np.arange(len(power_labels))  # the label locations
        width = 0.4  # the width of the bars

        ax1.bar(x - width / 2, sols_power_fo, width=width, color=color1)
        ax1b.bar(x + width / 2, sols_power_rt, width=width, color=color2)

        ax1.set_ylabel('Valid first-order solutions', color=color1)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1b.set_ylabel('Ray traceable solutions', color=color2)
        ax1b.tick_params(axis='y', labelcolor=color2)

        ax1.set_xticks(x)
        ax1.set_xticklabels(power_labels, rotation=45)

        # Orientation

        x = np.arange(len(orient_labels))  # the label locations
        width = 0.4  # the width of the bars

        ax2.bar(x - width / 2, sols_orient_fo, width=width, color=color1)
        ax2b.bar(x + width / 2, sols_orient_rt, width=width, color=color2)

        ax2.set_ylabel('Valid first-order solutions', color=color1)
        ax2.tick_params(axis='y', labelcolor=color1)
        ax2b.set_ylabel('Ray traceable solutions', color=color2)
        ax2b.tick_params(axis='y', labelcolor=color2)

        ax2.set_xticks(x)
        ax2.set_xticklabels(orient_labels, rotation=45)

        return

    def type_vs_spo(self):

        """


        Makes a violin plot of average spot size for different ray traceable
        solution spaces


        """

        # Initialize data structures

        spo = {}

        # Fill data structure with

        for sol in self.sols_rt:

            # Get power type

            power_str, orient_str = self.split_sol_type(sol.sol_type)
            if self.same_xy:
                power_str = self.abbrev_power_type(power_str, orient_str)

            # Fill dictionary with average spot size and average group EFL

            try:
                spo[power_str].append(sol.avg_spot_size)
            except KeyError:
                spo[power_str] = [sol.avg_spot_size]

        # Sort dictionary by median values, ascending

        num_type = len(spo)
        keys = np.empty(num_type, dtype='<U6')
        vals = np.empty(num_type)
        for i, (key, val) in enumerate(zip(spo.keys(), spo.values())):
            keys[i] = key
            vals[i] = np.median(val)

        order = np.argsort(vals)
        keys_sort = keys[order]

        # Initialize figure

        fig, ax = plt.subplots()

        # Plot violin plot of average spot size for different solution types

        vals = []
        labels = ['']
        median = np.empty(num_type)
        min_val = np.empty(num_type)
        max_val = np.empty(num_type)
        perc_25 = np.empty(num_type)
        perc_75 = np.empty(num_type)

        for i, key in enumerate(keys_sort):
            labels.append(key)
            vals.append(spo[key])
            min_val[i], \
            perc_25[i], \
            median[i], \
            perc_75[i],\
            max_val[i] = np.percentile(spo[key], [0, 25, 50, 75, 100])

        parts = ax.violinplot(vals, showextrema=False)

        inds = np.arange(num_type) + 1
        ax.vlines(inds, perc_25, perc_75, color='k', linestyle='-', lw=5)
        ax.vlines(inds, min_val, max_val, color='k', linestyle='-', lw=1)
        ax.scatter(inds, median, color='w', zorder=3)

        # Plot settings

        for pc in parts['bodies']:
            pc.set_facecolor('#37b34a')
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        ax.set_xticklabels(labels)
        ax.set_ylabel('Average spot size [mm]')
        ax.set_yscale('log')

        return

    def type_vs_packaging(self):

        """


        Makes violin plots of TTL, BFL, and clear aperture for different ray
        traceable solution spaces


        """

        # Initialize data structures

        ttl = {}
        bfl = {}
        ca = {}

        # Fill data structure with

        for sol in self.sols_rt:

            # Get power type

            power_str, orient_str = self.split_sol_type(sol.sol_type)
            if self.same_xy:
                power_str = self.abbrev_power_type(power_str, orient_str)

            # Fill dictionary with TTL values

            try:
                ttl[power_str].append(sol.ttl)
            except KeyError:
                ttl[power_str] = [sol.ttl]

            # Fill dictionary with BFL values

            try:
                bfl[power_str].append(sol.bfl)
            except KeyError:
                bfl[power_str] = [sol.bfl]

            # Fill dictionary with average clear aperture values

            try:
                ca[power_str].append(sol.avg_clear_aper)
            except KeyError:
                ca[power_str] = [sol.avg_clear_aper]

        # Initialize figure

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))
        plt.subplots_adjust(wspace=0.3)
        num_type = len(ttl)
        inds = np.arange(num_type) + 1

        # Plot TTL results

        # Sort dictionary by median values, ascending

        keys = np.empty(num_type, dtype='<U6')
        vals = np.empty(num_type)
        for i, (key, val) in enumerate(zip(ttl.keys(), ttl.values())):
            keys[i] = key
            vals[i] = np.median(val)

        order = np.argsort(vals)
        keys_sort = keys[order]

        # Plot violin plot of average spot size for different solution types

        vals = []
        labels = ['']
        median = np.empty(num_type)
        min_val = np.empty(num_type)
        max_val = np.empty(num_type)
        perc_25 = np.empty(num_type)
        perc_75 = np.empty(num_type)

        for i, key in enumerate(keys_sort):
            labels.append(key)
            vals.append(ttl[key])
            min_val[i], \
            perc_25[i], \
            median[i], \
            perc_75[i],\
            max_val[i] = np.percentile(ttl[key], [0, 25, 50, 75, 100])

        parts_ttl = ax1.violinplot(vals, showextrema=False)

        ax1.vlines(inds, perc_25, perc_75, color='k', linestyle='-', lw=5)
        ax1.vlines(inds, min_val, max_val, color='k', linestyle='-', lw=1)
        ax1.scatter(inds, median, color='w', zorder=3)

        # Plot settings

        for pc in parts_ttl['bodies']:
            pc.set_facecolor('#1b75bb')
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        ax1.set_xticklabels(labels)
        ax1.set_ylabel('Total track length [mm]')

        # Plot BFL results

        # Sort dictionary by median values, ascending

        keys = np.empty(num_type, dtype='<U6')
        vals = np.empty(num_type)
        for i, (key, val) in enumerate(zip(bfl.keys(), bfl.values())):
            keys[i] = key
            vals[i] = np.median(val)

        order = np.argsort(vals)
        keys_sort = keys[order]

        # Plot violin plot of average spot size for different solution types

        vals = []
        labels = ['']
        median = np.empty(num_type)
        min_val = np.empty(num_type)
        max_val = np.empty(num_type)
        perc_25 = np.empty(num_type)
        perc_75 = np.empty(num_type)

        for i, key in enumerate(keys_sort):
            labels.append(key)
            vals.append(bfl[key])
            min_val[i], \
            perc_25[i], \
            median[i], \
            perc_75[i], \
            max_val[i] = np.percentile(bfl[key], [0, 25, 50, 75, 100])

        parts_bfl = ax2.violinplot(vals, showextrema=False)

        ax2.vlines(inds, perc_25, perc_75, color='k', linestyle='-', lw=5)
        ax2.vlines(inds, min_val, max_val, color='k', linestyle='-', lw=1)
        ax2.scatter(inds, median, color='w', zorder=3)

        # Plot settings

        for pc in parts_bfl['bodies']:
            pc.set_facecolor('red')
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        ax2.set_xticklabels(labels)
        ax2.set_ylabel('Back focal length [mm]')

        # Plot CA results

        # Sort dictionary by median values, ascending

        keys = np.empty(num_type, dtype='<U6')
        vals = np.empty(num_type)
        for i, (key, val) in enumerate(zip(ca.keys(), ca.values())):
            keys[i] = key
            vals[i] = np.median(val)

        order = np.argsort(vals)
        keys_sort = keys[order]

        # Plot violin plot of average spot size for different solution types

        vals = []
        labels = ['']
        median = np.empty(num_type)
        min_val = np.empty(num_type)
        max_val = np.empty(num_type)
        perc_25 = np.empty(num_type)
        perc_75 = np.empty(num_type)

        for i, key in enumerate(keys_sort):
            labels.append(key)
            vals.append(ca[key])
            min_val[i], \
            perc_25[i], \
            median[i], \
            perc_75[i], \
            max_val[i] = np.percentile(ca[key], [0, 25, 50, 75, 100])

        parts_ca = ax3.violinplot(vals, showextrema=False)

        ax3.vlines(inds, perc_25, perc_75, color='k', linestyle='-', lw=5)
        ax3.vlines(inds, min_val, max_val, color='k', linestyle='-', lw=1)
        ax3.scatter(inds, median, color='w', zorder=3)

        # Plot settings

        for pc in parts_ca['bodies']:
            pc.set_facecolor('orange')
            pc.set_edgecolor('black')
            pc.set_alpha(1)

        ax3.set_xticklabels(labels)
        ax3.set_ylabel('Average element clear aperture [mm]')

        return

    def analyze_search_bound(self):

        """


        Makes histograms of the group EFL, TTL, and BFL Monte Carlo boundary
        conditions for both valid first-order and ray traceable solutions


        """

        # Initialize variables

        sols_list = [self.sols, self.sols_rt]
        num_sol_list = [self.num_sol, self.num_sol_rt]
        ylabel_str_list = ['Valid first-order ', 'Ray traceable ']

        num_bin = 100
        num_group = 5
        if self.var_type == 'CYL':
            num_group += 1

        buffer = 0.1
        buffer_1 = buffer * np.diff(self.config.efl_group_rng)
        buffer_2 = buffer * np.diff(self.config.ttl_rng)
        buffer_3 = buffer * np.diff(self.config.bfl_rng)

        # Initialize figure

        fig, axs = plt.subplots(len(sols_list), 3, figsize=(14, 7))
        plt.subplots_adjust(wspace=0.4, hspace=0.4)

        # Loop over solution sets

        for i, (sols, num_sol, ylabel_str) in enumerate(zip(sols_list,
                                                            num_sol_list,
                                                            ylabel_str_list)):

            # Initialize data structures

            group_efl = np.empty((num_sol, num_group))
            ttl = np.empty(num_sol)
            bfl = np.empty(num_sol)

            # Loop over all valid first-order solutions

            for j, sol in enumerate(sols):

                group_efl[j, :] = np.abs(sol.group_efl)
                ttl[j] = sol.ttl
                bfl[j] = sol.bfl

            # Plot histogram of group EFL, VL, and BFL for all valid first-order
            # solutions

            # Plot histograms

            axs[i, 0].hist(group_efl.flatten(), num_bin, color='#f8991d')
            axs[i, 1].hist(ttl, num_bin, color='#93268f')
            axs[i, 2].hist(bfl, num_bin, color='#129a48')

            # Plot settings

            axs[i, 0].set(xlabel='Group EFL, absolute value [mm]',
                    ylabel=ylabel_str + 'solution groups')
            axs[i, 1].set(xlabel='TTL [mm]', ylabel=ylabel_str + 'solutions')
            axs[i, 2].set(xlabel='BFL [mm]', ylabel=ylabel_str + 'solutions')

            axs[i, 0].set_xlim(self.config.efl_group_rng[0] - buffer_1,
                         self.config.efl_group_rng[1] + buffer_1)
            axs[i, 1].set_xlim(self.config.ttl_rng[0] - buffer_2,
                         self.config.ttl_rng[1] + buffer_2)
            axs[i, 2].set_xlim(self.config.bfl_rng[0] - buffer_3,
                         self.config.bfl_rng[1] + buffer_3)

            opacity = 0.25

            rect_1_lower = ptc.Rectangle((self.config.efl_group_rng[0], 0),
                                         -buffer_1, axs[i, 0].get_ylim()[1],
                                         facecolor='k', alpha=opacity)
            rect_1_upper = ptc.Rectangle((self.config.efl_group_rng[1], 0),
                                         buffer_1, axs[i, 0].get_ylim()[1],
                                         facecolor='k', alpha=opacity)
            rect_2_lower = ptc.Rectangle((self.config.ttl_rng[0], 0),
                                         -buffer_2, axs[i, 1].get_ylim()[1],
                                         facecolor='k', alpha=opacity)
            rect_2_upper = ptc.Rectangle((self.config.ttl_rng[1], 0),
                                         buffer_2, axs[i, 1].get_ylim()[1],
                                         facecolor='k', alpha=opacity)
            rect_3_lower = ptc.Rectangle((self.config.bfl_rng[0], 0),
                                         -buffer_3, axs[i, 2].get_ylim()[1],
                                         facecolor='k', alpha=opacity)
            rect_3_upper = ptc.Rectangle((self.config.bfl_rng[1], 0),
                                         buffer_3, axs[i, 2].get_ylim()[1],
                                         facecolor='k', alpha=opacity)

            axs[i, 0].add_patch(rect_1_lower)
            axs[i, 0].add_patch(rect_1_upper)
            axs[i, 1].add_patch(rect_2_lower)
            axs[i, 1].add_patch(rect_2_upper)
            axs[i, 2].add_patch(rect_3_lower)
            axs[i, 2].add_patch(rect_3_upper)

        return

    def group_efl_vs_spo(self):

        """


        Makes a scatter plot of group EFL vs average spot size for all ray
        traceable solution spaces


        """

        # Initialize data structures

        perf = {}

        # Fill data structure with

        for sol in self.sols_rt:

            # Get power type

            power_str, orient_str = self.split_sol_type(sol.sol_type)
            if self.same_xy:
                power_str = self.abbrev_power_type(power_str, orient_str)

            # Fill dictionary with average spot size and average group EFL

            try:
                perf[power_str][0].append(np.mean(np.abs(sol.group_efl)))
                perf[power_str][1].append(sol.avg_spot_size)
            except KeyError:
                perf[power_str] = [[np.mean(np.abs(sol.group_efl))],
                                   [sol.avg_spot_size]]

        # Initialize figure

        fig, ax = plt.subplots()

        # Plot scatter plot of average group EFL vs average spot size for
        # different solution types

        for key, vals in zip(perf.keys(), perf.values()):
            avg_group_efl = vals[0]
            avg_spot_size = vals[1]
            ax.scatter(avg_group_efl, avg_spot_size, s=4, label=key)

        # Plot settings

        ax.set(xlabel='Average group EFL, absolute value [mm]',
               ylabel='Average spot size [mm]')
        ax.set_yscale('log')
        ax.legend()

        return
