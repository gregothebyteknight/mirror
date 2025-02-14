import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from modules import analytical_pdf, marginal_r
from scipy.stats import beta, gamma
from const_case import sim, big_r_est
from gamma_case import sim_gamma, gamma_shift

sns.set_style('dark')


def pdf_const(big_r, n_sim):
    """
    Visualize the distribution of apparent radius {@r_set} in condition of a sphere of radius {@r_big}

    @r_big {int} - radius of a sphere
    @n_sim {int} - number of simulations
    """
    # SIMULATIONS
    r_set = sim(big_r, n_sim)

    # PLOT SIMULATIONS
    _, bin_edges, _ = plt.hist(r_set, density = True, alpha = 0.8,
    color = 'brown', label = "Simulation", bins = 100) # bins == size
    # Calculate midpoints of bin edges
    r_grid = (bin_edges[:-1] + bin_edges[1:]) / 2

    # ANALYTICAL SOLUTION 
    r_pdf =  analytical_pdf(r_set = r_grid, big_r = big_r)

    # BETA FIT
    a, b, loc, scale = beta.fit(data = r_set, floc = 0, fscale = big_r) # using the raw data
    r_beta = beta.pdf(x = r_grid, a = a, b = b, loc = loc, scale = scale)

    # COMPUTE MLE
    mle_beta = np.sum(beta.logpdf(x = r_set, a = a, b = b, 
                                  loc = loc, scale = scale)) # on the whole r dataset
    mle_analytical = np.sum(np.log(analytical_pdf(r_set = r_set, big_r = big_r))) # on the whole r dataset

    # PRINT ERRORS (AIC)
    print(f"AIC (Analytical): {round(-2 * mle_analytical, 2)}")
    print(f"AIC (Beta): {round(-2 * mle_beta + 2 * 2, 2)}")

    # PLOT CUSTOMIZATION
    plt.plot(r_grid, r_pdf, label = "Analytical PDF")
    plt.plot(r_grid, r_beta, label = "Beta Fit")

    plt.xlabel("Apparent Radius")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Distribution of Apparent Radius")
    plt.savefig("./images/pdf_const.png")

    return plt


def est_box_const(big_r, n_sim_set, n_rep):
    """
    Boxplot of estimated R for different number of simulations
    
    @big_r {int} - radius of a sphere
    @n_sim_set {list} - set of number of simulations
    @n_rep {int} - number of repetitions
    """
    # R ESTIMATION
    big_r_est_set = [
        [big_r_est(big_r, n_sim) for rep in range(n_rep)]
          for n_sim in n_sim_set
    ] # set of estimated R for different number of simulations

    _, axes = plt.subplots(1, 3, figsize = (15, 5), sharey = True)

    for ax, n_sim, estimates in zip(axes, n_sim_set, big_r_est_set):
        ax.boxplot(estimates,
                   boxprops = dict(color = "black"),  
                   whiskerprops = dict(color = "black"), 
                   capprops = dict(color = "black"), 
                   medianprops = dict(color = "green"))
        ax.axhline(y = big_r, color = 'red', linestyle = '--', 
                label = f"True value = {big_r}")
        ax.set_title(f'Number of Simulations = {n_sim}')
        ax.legend()

    axes[0].set_ylabel('Estimated R')
    plt.tight_layout()
    plt.savefig("./images/est_box_const.png")
    
    return plt


def pdf_gamma(theta, k, n_fol):
    """
    Visualize the distribution of apparent radius in condition of a sphere of radius R (gamma distribution)

    @theta {int} - shape of gamma distribution
    @k {int} - scale of gamma distribution
    """
    # SIMULATIONS
    r_set = sim_gamma(theta, k, n_fol)

    # ANALYTICAL PDF
    r_grid = np.linspace(0, np.max(r_set), 100) # x axis for plotting
    r_pdf = np.array([marginal_r(r, theta, k) for r in r_grid]) 

    # VISUALIZATION
    plt.hist(r_set, density = True, alpha = 0.8, label = "Simulation", bins = 100)
    plt.plot(r_grid, r_pdf, label = "Analytical PDF")

    # PLOT CUSTOMIZATION
    plt.xlabel("Apparent Radius")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Distribution of Apparent Radius")
    plt.savefig("./images/pdf_gamma.png")

    return plt


def viz_gamma_shift(theta, k, n_fol, stats = False):
    """
    Visualize the shift of gamma distribution

    @theta {int} - shape of gamma distribution
    @k {int} - scale of gamma distribution
    @n_fol {int} - number of follicles
    @stats {bool} - print the difference in parameters
    """
    # OBTAINING SHIFTED GAMMA PARAMETERS AND SIMULATIONS
    a, scale, r_set = gamma_shift(theta, k, n_fol)

    # GENERATION OF PDFs
    grid = np.linspace(0, np.max(r_set), 1000)

    r_big_gamma = gamma.pdf(grid, a = theta, scale = k) # using the grid
    r_gamma = gamma.pdf(grid, a = a, scale = scale) # reusing the grid

    if stats == True:
        # DIFFERENCE ANALYSIS
        print(f"Gamma of R: shape = {theta}, scale = {k}")
        print(f"Gamma of r: shape = {round(a, 2)}, scale = {round(scale, 2)}")
        print()
    
    # VISUALIZATION
    plt.plot(grid, r_big_gamma, label = "Gamma of R")
    plt.plot(grid, r_gamma, label = "Gamma of r")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Shift of Gamma Distribution")
    plt.savefig("./images/viz_gamma_shift.png")

    return plt


def viz_gamma_grid(theta_set, k_set, n_fol, mode = "Default"):
    """
    Visualize gamma distribution for different gamma parameters

    @theta_set {list} - set of theta values
    @k_set {list} - set of k values
    @n_fol {int} - number of follicles
    @mode {str} - mode of visualization
    """
    _, axes = plt.subplots(3, 3, figsize = (15, 15))
    for i, theta in enumerate(theta_set):
        for j, k in enumerate(k_set):
            ax = axes[i, j]
            plt.sca(ax) # set current axis

            if mode == "Shift": # Either default or shift
                print(f"Shift for θ={theta}, k={k}")
                viz_gamma_shift(theta, k, n_fol, stats = True) # still returnes plt
            else:
                pdf_gamma(theta, k, n_fol)

            ax.set_title(f"Initial Gamma Parameters: θ={theta}, k={k}")
            
    plt.tight_layout()
    plt.savefig("./images/viz_gamma_grid.png")
    
    return plt


def params_shift(params_set, mode = "theta"):
    """
    Visualize the initial parameters of gamma distribution 
    against the estimated ones
    
    @params_set {list} - set of parameters
    @mode {str} - mode of visualization
    """
    # Selection parameter for iteration
    mode_dict = {
        "theta": (
            lambda p: gamma_shift(theta = p, k = 2, n_fol = 1000),  # p is interpreted as θ
            "True θ",
            "Estimated θ",
            "Gamma Shift: θ Parameter"
        ),
        "k": (
            lambda p: gamma_shift(theta = 2, k = p, n_fol = 1000),  # p is interpreted as k
            "True k",
            "Estimated k",
            "Gamma Shift: k Parameter"
        )
    }

    func, xlabel, ylabel, title = mode_dict[mode]

    theta_set_prime = []
    k_set_prime = []
    for p in params_set:
        for i in range(30):
            res = func(p)
            theta_set_prime.append(res[0])
            k_set_prime.append(res[1])
    
    # Average the 30 repetitions for each parameter value.
    theta_set_prime = np.array(theta_set_prime).reshape(-1, 30).mean(axis=1)
    k_set_prime = np.array(k_set_prime).reshape(-1, 30).mean(axis=1)

    plt.scatter(params_set, theta_set_prime, label = "Estimated Value of θ")
    plt.scatter(params_set, k_set_prime, label = "Estimated Value of k")
    plt.plot(params_set, params_set, linestyle = '--', label = f"True Value of {mode}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.title(title)
    plt.savefig("./images/params_shift.png")

    return plt
