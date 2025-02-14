import numpy as np
from scipy.stats import gamma
from scipy.integrate import quad

def analytical_pdf(r_set, big_r): 
    """
    analytical probability density function for apparent radius (basically, likelihood - f(r|R))
    
    @r_set : set of r values
    @r_big : fixed R value
    """
    return r_set / (big_r * np.sqrt(big_r**2 - r_set**2) + 1e-10) # + 1e-10 for numerical stability


def marginal_r(r, theta, k):
    """
    Calculation of marginal distribution of r for a given R (basically, p(r))

    @r {int} - apparent radius
    @theta {int} - shape of gamma distribution
    @k {int} - scale of gamma distribution
    """
    integrand = lambda big_r: analytical_pdf(r, big_r) * gamma.pdf(big_r, a = theta, scale = k)
    integral, _ = quad(integrand, r, np.inf) # constraint: r <= R
    return integral



