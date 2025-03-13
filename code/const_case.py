import numpy as np


def sim(big_r, n_sim):
    """
    Simulation of {@n_sim} disections of r for sphere with a fixed R

    @big_r {int} - radius of a sphere
    @n_sim {int} - number of simulations per follicle
    """
    # SIMULATIONS
    d = np.random.uniform(low = 0, high = big_r, size = n_sim) 
    # where d is distance from center of a sphere to the section + symmetry assumption
    r_set = np.sqrt(big_r**2 - d**2) # apparent radius
    return r_set


def big_r_est(big_r, n_sim, stats = False):
    """
    Estimate R based on apparent radius

    @big_r {int} - real radius of a sphere
    @n_sim {int} - number of simulations per follicle
    """
    r_set = sim(big_r, n_sim)
    big_r_est = np.mean(r_set) * 4 / np.pi # mean estimation of R
    if stats == True:
        big_r_diff = np.abs(big_r_est - big_r) # difference between estimated and real R

        print(f"Estimated R: {round(big_r_est, 3)}")
        print(f"Difference: {round(big_r_diff, 3)}")
    return big_r_est
    