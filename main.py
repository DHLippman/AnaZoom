"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/20/20

"""

from designs import SystemConfig, AnamorphicZoom, Solutions
from monte_carlo import mc_search_cyl_var
from utilities import load_obj, save_obj, folder_exist
import numpy as np
import matplotlib.pyplot as plt
from time import time


def main():

    #sols = load_obj('1e8_cyl_var')
    #sols.demograph()
    #
    #return

    start = time()

    # Set system values

    config = SystemConfig(bfl_rng=np.array([35., 65.]),
                          ttl_rng=np.array([240., 365.]),
                          efl_group_rng=np.array([20., 500.]))

    # Perform Monte Carlo search for cylindrical variator solutions

    sols = mc_search_cyl_var(config, num_trial=1e8)
    print(sols)

    # print(sols)

    stop = time()

    print('Duration: {0:0.4f} hrs'.format((stop - start) / 3600))



    return


if __name__ == '__main__':
    main()
    plt.show()

