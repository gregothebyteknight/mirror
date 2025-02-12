import numpy as np
from scipy.stats import gamma

def sim_gamma(theta, k, n_fol):
    """ 
    Simulation of disections of r for each of {n_fol} follucules with radius R  
    distributed by a gamma law (R ~ Gamma(theta, k))

    @theta {int} - shape of gamma distribution
    @k {int} - scale of gamma distribution
    @n_fol {int} - number of follicles (for each of them we simulate one r)
    """
    r_set = [] # set of disections' r
    r_big_set = np.random.gamma(shape = theta, scale = k, 
    size = n_fol) # set of follicules' radii, or our prior p(R)
    
    # SIMULATIONS
    for fol in range(n_fol):
        r_big = r_big_set[fol]
        d = np.random.uniform(low = 0, high = r_big) # l from center of a sphere to the section
        r = np.sqrt(r_big**2 - d**2) # apparent radius, or raw r_small data
        r_set.append(r)

    return np.array(r_set)


def gamma_shift(theta, k, n_fol):
    """
    Suggesting simulated r as a shifted gamma distribution of R

    @theta {int} - shape of gamma distribution of R
    @k {int} - scale of gamma distribution of R
    @n_fol {int} - number of follicles (for each of them we simulate one r)
    """
    r_set = sim_gamma(theta, k, n_fol)

    # Gamma fit of r
    a, _, scale = gamma.fit(r_set, floc = 0) # fitting gamma distribution to r
    return a, scale, r_set


def big_r_gamma_est(theta, k, n_fol):
    r_set = sim_gamma(theta, k, n_fol)
    big_r_est = r_set * 4 / np.pi # mean estimation of R
    theta_est, _, k_est = gamma.fit(big_r_est) # fitting gamma distribution to R
    
    print(f"Estimated theta: {round(theta_est, 3)}")
    print(f"Estimated k: {round(k_est, 3)}")
    print(f"Theta difference: {round(theta_est - theta, 3)}")
    print(f"K difference: {round(k_est - k, 3)}")


