"""
Author:         David Henry Lippman
File:           main.py
Date created:   07/20/20
Date modified:  07/20/20

"""

from designs import SystemConfig
from monte_carlo import mc_search_cyl_var
import numpy as np
import matplotlib.pyplot as plt


def main():

    """

    # STEP 1: system configuration class

    # STEP 2: Monte Carlo search

    # STEP 3: Anamorphic zoom object

    # STEP 4: CODE V analysis

        # 4a) Create CODE V design

        # 4b) Check for ray trace failures at all zooms and fields

        # 4c) For successful ray traces, optimize and evaluate THO

    STEP 5: Demographics

        # 5a) Store demographics on solution types, number, ray traceability, and
            performance

        5b) Sankey plot

    STEP 6: Spherical variator

    """

    # Set system values

    config = SystemConfig(bfl_rng=np.array([35., 65.]),
                          ttl_rng=np.array([240., 365.]),
                          efl_group_rng=np.array([20., 150.]))

    # Perform Monte Carlo search for cylindrical variator solutions

    sols = mc_search_cyl_var(config, num_trial=16e6)

    print(sols)

    return


if __name__ == '__main__':
    main()
    plt.show()

