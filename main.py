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

    sols = load_obj('1e9_cyl_var_not_same_xy')

    for i, sol in enumerate(sols.sols_rt):
        sol.make_codev()
        ray_trace = sol.check_ray_trace()
        print(ray_trace)
        if not ray_trace:
            print(sol)

    return

    for sol in sols.sols_rt:
        print(sol.check_ray_trace())

    return

    # Set system values

    config = SystemConfig(bfl_rng=np.array([35., 65.]),
                          ttl_rng=np.array([240., 360.]),
                          efl_group_rng=np.array([20., 500.]))

    # Perform Monte Carlo search for cylindrical variator solutions

    sols = mc_search_cyl_var(config, num_trial=1e9, same_xy=False)
    
    print(sols)

    return





if __name__ == '__main__':
    main()
    plt.show()

