"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/28/20

"""

from designs import SystemConfig, AnamorphicZoom, Solutions
from monte_carlo import mc_search_cyl_var, mc_search_sph_var
from utilities import load_obj, save_obj, folder_exist
import numpy as np
import matplotlib.pyplot as plt


def main():

    # Set system values

    config = SystemConfig(bfl_rng=np.array([35., 65.]),
                          ttl_rng=np.array([240., 365.]),
                          efl_group_rng=np.array([20., 500.]))

    # Perform Monte Carlo search for solutions

    # sols = mc_search_cyl_var(config, num_trial=1e6, same_xy=False)
    sols = mc_search_sph_var(config, num_trial=1e6)

    print(sols)

    return


if __name__ == '__main__':
    main()
    plt.show()

