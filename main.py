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

    sys = SystemSetup()

    # TODO: upload to Github

    return


class SystemSetup:

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

        ttl_rng:            (2,) array of total track length range mm]; defined
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

        # TEMP: plot

        fig, axs = plt.subplots(1, self.num_zoom)

        for z in range(self.num_zoom):
            axs[z].scatter(np.degrees(self.xan[:, z]),
                           np.degrees(self.yan[:, z]))
            axs[z].set(xlabel='X field angle [deg]',
                       ylabel='Y field angle [deg]',
                       title='EFX = {0:0.2f} mm\nEFY = {1:0.2f} mm'
                             .format(self.efx[z], self.efy[z]))
            axs[z].set_aspect('equal')


        return


if __name__ == '__main__':
    main()
    plt.show()

