"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/20/20

"""

from designs import SystemConfig
from monte_carlo import mc_search_cyl_var
from utilities import load_obj
import numpy as np
import matplotlib.pyplot as plt


def main():

    # Set system values

    config = SystemConfig(bfl_rng=np.array([35., 65.]),
                          ttl_rng=np.array([240., 365.]),
                          efl_group_rng=np.array([20., 500.]))

    # Perform Monte Carlo search for cylindrical variator solutions

    sols = mc_search_cyl_var(config, num_trial=1e8)

    print(sols)

    return


if __name__ == '__main__':
    main()
    plt.show()

