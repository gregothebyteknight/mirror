# Simulation of Spheres and Tissue Structures

This project encompasses a series of simulations and analytical derivations related to the geometric properties of spheres and the generation of 3D biological tissue structures. The primary goals are:
1.  To understand the distribution of apparent 2D cross-section radii ($r$) obtained from slicing spheres of true 3D radii ($R$).
2.  To infer the true sphere radii ($R$), or parameters of their distribution, from observed cross-section radii $r$.
3.  To simulate realistic 3D structures of nodules and vascular networks using various stochastic models.

The repository includes Python code for sphere-related simulations and R code for tissue generation, along with detailed notes and visualization scripts.

## Theoretical Background

This project simulates and analyzes the geometric properties of spheres, focusing on inferring true 3D sphere radii ($R$) from observed 2D cross-section radii ($r$). It also includes modules for simulating 3D biological structures like nodules and vessels.

### Single Sphere Case

Consider a single sphere of a fixed, known radius $R$. When this sphere is randomly sliced by a plane, a circular cross-section (dissection) is formed with a radius $r \le R$.

The key assumptions are:
1.  The distance $d$ from the center of the sphere to the slicing plane is uniformly distributed over the interval $[0, R]$. We consider only one hemisphere due to symmetry, so $d \sim U(0, R)$. The probability density function (PDF) for $d$ is $f_D(d) = \frac{1}{R}$ for $0 \le d \le R$.
2.  The radius of the dissection $r$ is related to $R$ and $d$ by the Pythagorean theorem: $r^2 + d^2 = R^2$. Thus, $r = \sqrt{R^2 - d^2}$. Correspondingly, $d = \sqrt{R^2 - r^2}$.

To find the distribution of the dissection radius $r$, we use the method of transformations.
The PDF of $r$ given $R$, denoted $f(r|R)$, can be derived as:
$$f(r|R) = f_D(d) \left| \frac{\partial d}{\partial r} \right|$$
Since $d = \sqrt{R^2 - r^2}$, the derivative is:
$$\frac{\partial d}{\partial r} = \frac{-r}{\sqrt{R^2 - r^2}}$$
And its absolute value is $\left| \frac{-r}{\sqrt{R^2 - r^2}} \right| = \frac{r}{\sqrt{R^2 - r^2}}$ (since $r \ge 0$).

Therefore, the analytical PDF for the dissection radius $r$ (for $0 \le r \le R$) is:
$$f(r|R) = \frac{1}{R} \cdot \frac{r}{\sqrt{R^2 - r^2}}$$
This formula is implemented in `modules.py` as `analytical_pdf(r_set, big_r)`.

### Set of Spheres (Gamma Distributed Radii)

In a more complex scenario, we consider a population of spheres where the true radii $R$ are not fixed but are themselves drawn from a probability distribution. This project assumes that the true radii $R$ follow a Gamma distribution: $R \sim \text{Gamma}(\alpha, \beta)$, where $\alpha$ is the shape parameter and $\beta$ is the scale parameter. The PDF of $R$ is denoted $f_R(R; \alpha, \beta)$.

When each of these spheres is dissected, we observe a set of dissection radii $r$. The distribution of these observed $r$ values is a marginal distribution, obtained by integrating the conditional distribution $f(r|R)$ over all possible values of $R$, weighted by the probability of each $R$:
$$f(r; \alpha, \beta) = \int_{R=r}^{\infty} f(r|R) f_R(R; \alpha, \beta) dR$$
The lower limit of integration is $r$ because a dissection radius $r$ cannot be larger than the true sphere radius $R$.
This integral gives the marginal probability density of observing a dissection radius $r$, given the parameters $\alpha, \beta$ of the Gamma distribution for $R$. This is implemented in `modules.py` as `marginal_r(r, theta, k)` (where `theta` corresponds to $\alpha$ and `k` to $\beta$).

### R Inference

A primary goal of the project is to estimate the true sphere radius $R$ (or its distribution parameters) from the observed dissection radii $r$.

#### Single R Case (Constant R)

If we have a set of observed dissection radii $r_1, r_2, \dots, r_n$ from a single sphere of unknown (but constant) radius $R$, we can estimate $R$.
The expected value of $r$, $\mathbb{E}(r)$, can be calculated using $f(r|R)$:
$$\mathbb{E}(r) = \int_0^R r \cdot f(r|R) dr = \int_0^R r \cdot \frac{1}{R} \frac{r}{\sqrt{R^2 - r^2}} dr$$
This integral evaluates to:
$$\mathbb{E}(r) = \frac{1}{R} \left[ \frac{R^2}{2} \arcsin\left(\frac{r}{R}\right) - \frac{r}{2}\sqrt{R^2-r^2} \right]_0^R = \frac{1}{R} \left( \frac{R^2}{2} \cdot \frac{\pi}{2} \right) = \frac{\pi R}{4}$$
From this, we can derive an estimator for $R$ based on the sample mean of $r$ (denoted $\bar{r}$):
$$\hat{R} = \bar{r} \cdot \frac{4}{\pi}$$
This estimator is implemented in `const_case.py` within the `big_r_est` function.

#### R Distribution Case (Gamma Distributed R)

When the true radii $R$ follow a Gamma distribution $R \sim \text{Gamma}(\alpha, \beta)$, and we only observe one dissection $r_i$ from each sphere $R_i$, estimating the parameters $\alpha$ and $\beta$ of the $R$ distribution is challenging.
A direct application of the $\hat{R}_i = r_i \cdot \frac{4}{\pi}$ formula for each $r_i$ would yield highly noisy estimates for individual $R_i$'s, as each $r_i$ is a single sample from its conditional distribution $f(r|R_i)$.

The project explores the relationship between the parameters of the Gamma distribution fitted to the observed $r$ values (let's say $r \sim \text{Gamma}(\alpha', \beta')$) and the original parameters $\alpha, \beta$ of the $R$ distribution. The idea is to find functions $g_1, g_2$ such that $\alpha = g_1(\alpha')$ and $\beta = g_2(\beta')$. This is investigated empirically by simulating $R$ values from a known Gamma distribution, generating corresponding $r$ values, fitting a Gamma distribution to these $r$ values, and then trying to model the transformation between the true and estimated parameters. This is explored in `gamma_case.py` and visualized in `viz.py` (e.g., `params_shift` function).

### Simulation of Tissues

Beyond spherical dissections, the project includes modules for simulating 3D biological tissue structures.

#### Nodules

Nodule simulation aims to generate clusters of points in 3D space. This is achieved using a **Poisson cluster process**:
1.  **Cluster Centroids**: The locations of nodule centers are first generated from a 3D Poisson point process within a defined volume.
2.  **Points around Centroids**: For each centroid, a number of points (cells) are generated around it. The distribution of these points relative to the centroid can be controlled, for example, by drawing a radius for each cluster (e.g., from an Exponential distribution) and then distributing points within that sphere according to another Poisson process.
This is implemented in `code/nodule.r`.

#### Vessels

Vessel simulation aims to create branching, tree-like structures in 3D. Two main approaches are mentioned:
1.  **Brownian Random Walk**: A simpler model where vessel paths are generated using a 3D Brownian motion, with branching occurring probabilistically. This is explored in `code/vessel.py`.
2.  **Kent Distribution Model**: A more sophisticated approach for directional control of branching. The Kent distribution is a directional distribution on the unit sphere, allowing for more realistic branching patterns.
    *   An iterative process generates a tree of vessel segment centroids.
    *   At each step, a new segment direction can be drawn from a Kent distribution.
    *   Branching can occur with a certain probability at the tip of existing segments.
    *   Once the vessel skeleton (centroids) is formed, points representing vessel cells are distributed around these segments, often using a Poisson process within a certain radius of the segment.
This model is primarily implemented in `code/vessel.r`, utilizing helper functions from `code/module.r`.

---

## Code Structure and Usage

The project is organized into Python scripts, R scripts, a Jupyter Notebook for orchestrating Python simulations, and data/note files.

### Python Modules (located in `code/`)

The Python part of the project focuses on the simulation of sphere dissections and inference of sphere radii.

**`const_case.py`**
*   **Purpose**: Simulates dissections for a single sphere with a constant radius $R$ and provides functions to estimate $R$.
*   **Key Functions**:
    *   `sim(big_r, n_sim)`: Simulates `n_sim` dissection radii $r$ for a sphere of radius `big_r` ($R$). Returns a NumPy array of $r$ values.
    *   `big_r_est(big_r, n_sim, stats=False)`: Estimates the true radius $R$ using the mean estimator formula ($\hat{R} = \text{mean}(r) \cdot 4/\pi$) based on `n_sim` simulated dissections. If `stats=True`, prints the estimate and error.
*   **Usage**:
    ```python
    from const_case import sim, big_r_est

    true_R = 10.0
    num_simulations = 1000

    # Simulate dissection radii
    r_values = sim(true_R, num_simulations)

    # Estimate R
    estimated_R = big_r_est(true_R, num_simulations, stats=True)
    print(f"Simulated r values (first 5): {r_values[:5]}")
    print(f"Estimated R: {estimated_R}")
    ```

**`gamma_case.py`**
*   **Purpose**: Simulates dissections when true sphere radii $R$ follow a Gamma distribution ($R \sim \text{Gamma}(\theta, k)$). Includes functions for fitting Gamma distributions to observed $r$ values and attempting to estimate the parameters of the original $R$ distribution.
*   **Key Functions**:
    *   `sim_gamma(theta, k, n_fol)`: Simulates `n_fol` dissection radii $r$, where each corresponding true radius $R$ is drawn from $\text{Gamma}(\text{theta}, \text{k})$. Returns a NumPy array of $r$ values.
    *   `gamma_shift(theta, k, n_fol)`: Simulates $r$ values using `sim_gamma`, then fits a Gamma distribution to these $r$ values. Returns the fitted shape (`a`), scale (`scale`), and the raw $r$ values.
    *   `big_r_gamma_est(theta, k, n_fol)`: Simulates $r$ values, applies the mean estimator ($R_{est} = r \cdot 4/\pi$) to each $r$ to get estimates of $R$, then fits a Gamma distribution to these $R_{est}$ values. Prints estimated vs. true parameters.
*   **Usage**:
    ```python
    from gamma_case import sim_gamma, gamma_shift, big_r_gamma_est

    shape_R = 2.0
    scale_R = 3.0
    num_follicles = 1000

    # Simulate r values where R ~ Gamma(shape_R, scale_R)
    r_values_gamma = sim_gamma(shape_R, scale_R, num_follicles)

    # Fit Gamma to r values
    fitted_shape_r, fitted_scale_r, _ = gamma_shift(shape_R, scale_R, num_follicles)
    print(f"Fitted Gamma params for r: shape={fitted_shape_r:.2f}, scale={fitted_scale_r:.2f}")

    # Estimate parameters of R distribution
    big_r_gamma_est(shape_R, scale_R, num_follicles)
    ```

**`modules.py`**
*   **Purpose**: Contains core mathematical/analytical functions used by other Python modules.
*   **Key Functions**:
    *   `analytical_pdf(r_set, big_r)`: Calculates the analytical PDF $f(r|R) = \frac{1}{R} \cdot \frac{r}{\sqrt{R^2 - r^2}}$ for a given set of $r$ values (`r_set`) and a fixed true radius `big_r` ($R$).
    *   `marginal_r(r, theta, k)`: Calculates the marginal PDF $f(r)$ for a dissection radius $r$, given that the true radii $R$ follow $\text{Gamma}(\text{theta}, \text{k})$. This involves numerical integration.
*   **Usage**: These functions are typically called internally by `viz.py` or other simulation modules.

**`viz.py`**
*   **Purpose**: Handles visualization of simulation results from `const_case.py` and `gamma_case.py`. Generates and saves various plots.
*   **Key Functions**:
    *   `pdf_const(big_r, n_sim)`: Plots the histogram of simulated $r$ values (from `const_case.sim`), overlays the analytical PDF and a fitted Beta distribution. Prints AIC values for comparison. Saves plot to `images/pdf_const.pdf`.
    *   `est_box_const(big_r, n_sim_set, n_rep)`: Generates boxplots comparing $R$ estimates for different numbers of simulations (`n_sim_set`), repeated `n_rep` times. Saves plot to `images/est_box_const.pdf`.
    *   `pdf_gamma(theta, k, n_fol)`: Plots the histogram of simulated $r$ values (from `gamma_case.sim_gamma`) and overlays the analytical marginal PDF. Saves plot to `images/pdf_gamma.pdf`.
    *   `viz_gamma_shift(theta, k, n_fol, stats=False)`: Plots the PDF of the true $R \sim \text{Gamma}(\text{theta}, \text{k})$ against the PDF of the Gamma distribution fitted to the simulated $r$ values. Saves plot to `images/viz_gamma_shift.pdf`.
    *   `viz_gamma_grid(theta_set, k_set, n_fol, mode)`: Creates a grid of plots, calling either `pdf_gamma` or `viz_gamma_shift` for different combinations of `theta` and `k`. Saves plot to `images/viz_gamma_grid.pdf`.
    *   `params_shift(params_set, mode)`: Visualizes the relationship between true Gamma parameters ($\theta$ or $k$) and their estimated values derived from fitting the $r$ distribution. Saves plot to `images/params_shift.pdf`.
*   **Usage**: Typically run via `main.ipynb`, which calls these functions with appropriate parameters.

**`vessel.py`**
*   **Purpose**: Implements a 3D branching Brownian motion model for simulating vessel-like structures. This seems to be one of the explored models for vessel generation.
*   **Key Functions**:
    *   `walk()`: Simulates and plots a single 1D Brownian motion path with smoothing.
    *   `tree_walk()`: Simulates a 1D branching Brownian motion, where paths can split. Plots multiple smoothed branches.
    *   `space_walk()`: Simulates a 3D branching Brownian motion. Plots smoothed 3D paths.
*   **Usage**: Can be run as a script (`if __name__ == "__main__":`).
    ```bash
    python code/vessel.py
    ```
    This will generate and display plots.

**`main.ipynb`**
*   **Purpose**: A Jupyter Notebook that serves as the primary interface for running the Python-based sphere dissection simulations and visualizations.
*   **Structure**:
    *   Imports necessary functions from `const_case.py`, `gamma_case.py`, and `viz.py`.
    *   Contains cells for setting parameters (e.g., `big_r`, `n_sim`, `theta`, `k`, `n_fol`, sets for grid plots).
    *   Calls visualization functions from `viz.py` to generate and save plots.
*   **Usage**:
    1.  Open `code/main.ipynb` in a Jupyter Notebook environment.
    2.  Modify parameters in the code cells as needed.
    3.  Run cells sequentially to execute simulations and generate visualizations. Plots will be saved to the `images/` directory.

### R Scripts (located in `code/` and `archive/`)

The R scripts focus on the simulation of 3D tissue structures (nodules and vessels).

**`code/nodule.r`**
*   **Purpose**: Generates 3D "nodule" point patterns using a Poisson cluster process.
*   **Main Function**: `nodule(int, ord, viz = FALSE, save = FALSE)`
    *   `int`: Overall point intensity.
    *   `ord`: Simulation index (for labeling).
    *   `viz`: If `TRUE`, renders a 3D WebGL HTML visualization (`images/nodule.html`).
    *   `save`: If `TRUE`, writes centered coordinates to `data/nodule_loc.csv`.
*   **Usage**:
    ```R
    source("code/nodule.r")
    # Example: Generate nodules, visualize, and save
    nodule_data <- nodule(int = 7e4, ord = 1, viz = TRUE, save = TRUE)
    ```

**`code/vessel.r`**
*   **Purpose**: Simulates 3D "vessel" point patterns, likely using a branching model based on the Kent distribution for trajectory growth and Poisson-distributed points for vessel walls.
*   **Main Function**: `vessel(int, ord, viz = FALSE, save = FALSE)`
    *   `int`: Poisson intensity for wall-points.
    *   `ord`: Simulation identifier.
    *   `viz`: If `TRUE`, renders a 3D WebGL HTML visualization (`images/vessel_combined.html`).
    *   `save`: If `TRUE`, writes centered coordinates to `data/vessel_loc.csv`.
*   **Dependencies**: Uses helper functions from `code/module.r` (e.g., for Kent distribution points, spatial Poisson processes).
*   **Usage**:
    ```R
    source("code/vessel.r")
    # Example: Generate vessels, visualize
    vessel_data <- vessel(int = 5e-3, ord = 1, viz = TRUE, save = FALSE)
    ```

**`code/module.r`**
*   **Purpose**: Contains shared helper functions used by other R scripts, particularly `vessel.r`.
*   **Key Functions (inferred from notes)**:
    *   `kent_point()`: Generates points based on the Kent distribution for directional vessel growth.
    *   `pois_in_space()`: Simulates points in a 3D space (e.g., parallelepiped) according to a Poisson process, used for populating vessel walls.
    *   `rotation()`: Performs 3D rotations, likely used to align vessel segments.
*   **Usage**: Sourced by other R scripts. Not typically run directly.

**`code/cultivate.r`**
*   **Purpose**: Appears to be a master R script for running multiple tissue simulations (either "nodule" or "vessel" type), saving their outputs, and potentially generating summary plots (e.g., boxplot of cell counts).
*   **Functionality (inferred from `notes/simulations.md`)**:
    *   Sets parameters like intensity, number of simulations, and simulation type ("nodule" or "vessel").
    *   Loops through replicate simulations, calling `nodule.r` or `vessel.r`.
    *   Saves individual simulation outputs (e.g., coordinates to `data/sim_<i>/selected_cell_coordinates.csv`).
    *   May generate summary plots (e.g., `images/<type>_cell_counts_boxplot.png`).
*   **Usage**:
    ```bash
    Rscript code/cultivate.r
    ```
    Parameters like simulation type and number of runs are typically set inside the script.

**`archive/test.r`**
*   **Purpose**: Likely an older R script used for testing specific functionalities or experimental code. Its exact current relevance might be limited.

### Data Files (located in `data/`)

*   **`data/nodule_loc.csv`**: Stores X, Y, Z coordinates and cluster/cell labels for points generated by `code/nodule.r` when `save=TRUE`.
*   **`data/vessel_loc.csv`**: Stores X, Y, Z coordinates and cluster/cell labels for points generated by `code/vessel.r` when `save=TRUE`.
*   The `code/cultivate.r` script might also generate simulation-specific CSV files in subdirectories like `data/sim_<i>/`.

---

## Setup and Installation

### Python Environment

The Python scripts rely on several common scientific computing libraries. It's recommended to use a virtual environment (e.g., `venv` or `conda`) for managing dependencies.

**Dependencies:**
*   Python (3.6+ recommended)
*   NumPy: For numerical operations.
*   SciPy: For scientific and technical computing (used for `gamma.fit`, `quad` integration, `savgol_filter`).
*   Matplotlib: For plotting.
*   Seaborn: For enhanced statistical visualizations.
*   Jupyter Notebook or JupyterLab: For running `main.ipynb`.

**Installation:**

1.  **Create and activate a virtual environment (optional but recommended):**
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```
2.  **Install Python packages:**
    You can install the packages using pip:
    ```bash
    pip install numpy scipy matplotlib seaborn jupyterlab
    ```
    (If you prefer Jupyter Notebook over JupyterLab, replace `jupyterlab` with `notebook`).

### R Environment

The R scripts require R to be installed, along with several packages from CRAN.

**Dependencies:**
*   R (latest version recommended)
*   `rgl`: For 3D interactive visualizations (nodules, vessels).
*   `htmlwidgets`: For saving `rgl` visualizations as HTML files.
*   `fifer`: (Mentioned in notes for `vessel.r`, though its specific use isn't detailed in the provided snippets. May be for specific statistical functions or data manipulation).
*   `spatstat`: (Mentioned in notes for `vessel.r`, likely for spatial statistics or point process generation).
*   `compositions`: (While not explicitly in `simulations.md`, the `archive/test.r` uses `library(compositions)`, so it might be relevant for some parts or older code).

**Installation:**

1.  **Install R:** Download and install R from [CRAN (The Comprehensive R Archive Network)](https://cran.r-project.org/).
2.  **Install R packages:**
    Open an R console and run the following commands:
    ```R
    install.packages(c("rgl", "htmlwidgets", "fifer", "spatstat", "compositions"))
    ```

---

## Running Simulations and Visualizations

Parameters for simulations (e.g., sphere radius, number of simulations, Gamma distribution parameters, intensities for tissue models) are generally hard-coded within the respective scripts or the `main.ipynb` notebook. You will need to edit these files directly to change simulation conditions.

### Python Simulations (Sphere Dissections)

The primary way to run the Python-based sphere dissection simulations and generate visualizations is through the Jupyter Notebook:

1.  **Start JupyterLab or Jupyter Notebook:**
    If you installed JupyterLab:
    ```bash
    jupyter lab
    ```
    Or for Jupyter Notebook:
    ```bash
    jupyter notebook
    ```
    This will open a new tab in your web browser.
2.  **Navigate and open `code/main.ipynb`**.
3.  **Modify Parameters:** Inside the notebook, locate the cells where parameters are defined (e.g., `big_r`, `n_sim` for the constant case; `theta`, `k`, `n_fol` for the Gamma case; `n_sim_set`, `theta_set`, `k_set` for batch visualizations). Adjust these values as needed.
4.  **Run Cells:** Execute the notebook cells sequentially (e.g., using "Run" -> "Run All Cells" or by running cells individually with Shift+Enter).
    *   Simulations will be performed.
    *   Visualization functions from `viz.py` will be called.
    *   Output plots will be saved to the `images/` directory (e.g., `pdf_const.pdf`, `est_box_const.pdf`, `pdf_gamma.pdf`, etc.).
    *   Some functions like `big_r_est` (when `stats=True`) or `big_r_gamma_est` will print results directly in the notebook.

The `code/vessel.py` script (Brownian motion vessel model) can be run directly as a Python script:
```bash
python code/vessel.py
```
This will typically display Matplotlib plots on screen.

### R Simulations (Tissue Structures)

The R scripts for simulating nodules and vessels can be run from the command line using `Rscript` or interactively within an R environment.

**1. Individual Nodule/Vessel Simulations:**

You can source `nodule.r` or `vessel.r` in an R console and call their main functions:

```R
# In an R console
# For Nodules:
source("code/nodule.r")
nodule_data <- nodule(int = 7e4, ord = 1, viz = TRUE, save = TRUE)
# This will generate images/nodule.html and data/nodule_loc.csv

# For Vessels:
source("code/vessel.r") # This will also source module.r if needed
vessel_data <- vessel(int = 5e-3, ord = 1, viz = TRUE, save = FALSE)
# This will generate images/vessel_combined.html
```
Modify the parameters (`int`, `ord`, `viz`, `save`) as needed.

**2. Batch Tissue Simulations (`cultivate.r`):**

The `code/cultivate.r` script is designed to run multiple simulations and save their outputs.
*   **Modify Parameters within `cultivate.r`**: Open `code/cultivate.r` and adjust variables like `type` ("nodule" or "vessel"), `int_nod` or `int_vess` (intensities), and `num_sim` (number of replicates).
*   **Run from Command Line**:
    ```bash
    Rscript code/cultivate.r
    ```
    This will:
    *   Create output folders (e.g., `data/sim_<i>/`).
    *   Run the chosen simulation type `num_sim` times.
    *   Save coordinate data for each simulation.
    *   Potentially generate a summary boxplot (e.g., `images/<type>_cell_counts_boxplot.png`).

**Important Notes for R Scripts:**
*   Ensure that the working directory is set to the root of the project so that relative paths like `code/module.r` or `../data/` resolve correctly. If running from an R console, you might need to use `setwd()` to the project root. When using `Rscript` from the project root, paths should generally work as written in the scripts.
*   For visualizations (`viz = TRUE`), R will attempt to open an `rgl` window and save an HTML widget. This might require a graphical environment if not running on a local machine.

---

## Examples

Here are a couple of quick examples to get started.

### Example 1: Python - Single Sphere Simulation and R Estimation

This example uses Python to simulate dissections from a single sphere and estimate its radius. It mirrors parts of what `main.ipynb` does.

1.  Ensure your Python environment is set up and you are in the project's root directory.
2.  Create a Python script, e.g., `run_single_sphere_example.py`:
    ```python
    from code.const_case import sim, big_r_est
    from code.viz import pdf_const, est_box_const # For visualization
    import matplotlib.pyplot as plt # Required by viz module for showing plots if not in notebook

    # Parameters
    true_R = 7.5
    num_simulations_for_pdf = 2000

    # 1. Simulate and visualize PDF of apparent radii
    print(f"Running PDF visualization for R = {true_R} with {num_simulations_for_pdf} simulations...")
    plt.figure() # Create a new figure for the plot
    pdf_const(big_r=true_R, n_sim=num_simulations_for_pdf)
    plt.suptitle(f"PDF for R={true_R}, n_sim={num_simulations_for_pdf}") # Add a main title
    plt.show() # Display the plot
    print(f"Plot saved to images/pdf_const.pdf (or similar name from viz.py)")

    # 2. Estimate R and check convergence with boxplots
    print("\nRunning R estimation boxplots...")
    n_sim_set_for_boxplot = [50, 200, 1000] # Number of simulations per estimate
    num_repetitions_for_boxplot = 50       # Number of times R is estimated for each n_sim in set

    plt.figure() # Create a new figure for the plot
    est_box_const(big_r=true_R, n_sim_set=n_sim_set_for_boxplot, n_rep=num_repetitions_for_boxplot)
    plt.suptitle(f"R Estimation Boxplots for True R={true_R}") # Add a main title
    plt.show() # Display the plot
    print(f"Plot saved to images/est_box_const.pdf (or similar name from viz.py)")

    # 3. Single R estimation
    num_simulations_for_est = 500
    estimated_R = big_r_est(true_R, num_simulations_for_est, stats=True)
    print(f"\nFor a single run with {num_simulations_for_est} simulations:")
    print(f"True R: {true_R}, Estimated R: {estimated_R:.3f}")
    ```
3.  Run the script from your terminal:
    ```bash
    python run_single_sphere_example.py
    ```
    This will:
    *   Print AIC values and save `images/pdf_const.pdf`.
    *   Save `images/est_box_const.pdf`.
    *   Print the single R estimation results.
    *   Display the plots on screen.

### Example 2: R - Nodule Simulation

This example uses R to generate a 3D nodule structure, visualize it, and save the coordinates.

1.  Ensure your R environment is set up and you are in the project's root directory.
2.  Run the following in an R console or save it as an R script (e.g., `run_nodule_example.R`) and run `Rscript run_nodule_example.R`:
    ```R
    # Set working directory to project root if not already there
    # setwd("/path/to/your/SimulationOfSpheres/")

    source("code/nodule.r")

    # Parameters
    intensity <- 5e4 # Lower intensity for a quicker example
    simulation_id <- "example1"

    cat("Simulating nodules...\n")
    nodule_output <- nodule(
      int = intensity,
      ord = simulation_id,
      viz = TRUE,  # Generate images/nodule_example1.html (name might vary based on actual nodule.r implementation)
      save = TRUE  # Save data/nodule_loc_example1.csv (name might vary)
    )

    cat("Nodule simulation complete.\n")
    if (TRUE) { # Assuming viz = TRUE
      cat("Visualization saved to an HTML file in the 'images' directory (e.g., images/nodule.html or similar).\n")
    }
    if (TRUE) { # Assuming save = TRUE
      cat("Nodule coordinates saved to a CSV file in the 'data' directory (e.g., data/nodule_loc.csv or similar).\n")
    }

    # To inspect the returned data (if not saving, or to see it anyway)
    # print(head(nodule_output))
    ```
    This will:
    *   Simulate nodule points.
    *   Create an interactive 3D HTML visualization in the `images` directory.
    *   Save the nodule coordinates to a CSV file in the `data` directory.
    *   Print messages to the console.

---

## Directory Structure

```
SimulationOfSpheres/
├── .gitignore
├── LICENSE
├── README.md                 # This file
├── archive/
│   └── test.r                # Archived R test script
├── code/
│   ├── __pycache__/          # Python cache files
│   ├── const_case.py         # Python: Single sphere (constant R) simulation
│   ├── cultivate.r           # R: Master script for batch tissue simulations
│   ├── gamma_case.py         # Python: Sphere set (Gamma R) simulation
│   ├── main.ipynb            # Jupyter Notebook: Main runner for Python simulations
│   ├── module.r              # R: Helper modules for R scripts (e.g., Kent distribution)
│   ├── modules.py            # Python: Core mathematical/analytical functions
│   ├── nodule.r              # R: Nodule simulation
│   ├── vessel.py             # Python: Vessel simulation (Brownian motion model)
│   ├── vessel.r              # R: Vessel simulation (Kent distribution model)
│   └── viz.py                # Python: Visualization scripts
├── data/
│   ├── nodule_loc.csv        # Example output: Nodule coordinates
│   └── vessel_loc.csv        # Example output: Vessel coordinates
│   └── sim_<i>/              # Potential subdirectories from cultivate.r
│       └── selected_cell_coordinates.csv
├── images/                   # Output directory for plots and visualizations
│   ├── pdf_const.pdf         # Output from viz.py: PDF for constant R case
│   ├── est_box_const.pdf     # Output from viz.py: Boxplot of R estimates
│   ├── pdf_gamma.pdf         # Output from viz.py: PDF for Gamma R case
│   ├── viz_gamma_shift.pdf   # Output from viz.py: Gamma distribution shift plot
│   ├── viz_gamma_grid.pdf    # Output from viz.py: Grid of Gamma visualizations
│   ├── params_shift.pdf      # Output from viz.py: Gamma parameter shift plot (Note: existing images/params_shift.png might be an older version)
│   ├── nodule.html           # Output from nodule.r: 3D Nodule visualization
│   ├── vessel_combined.html  # Output from vessel.r: 3D Vessel visualization
│   ├── *.png                 # Other existing images (e.g., est_box_const.png) or R-generated images
├── logs/
│   └── logs_slice_2025-03-19.txt # Example log file
├── notes/
│   ├── geometry.md           # Detailed notes on sphere geometry simulations
│   └── simulations.md        # Detailed notes on tissue simulations (R scripts)
└── run_single_sphere_example.py # Example script created for README
└── run_nodule_example.R         # Example script created for README
```

---
