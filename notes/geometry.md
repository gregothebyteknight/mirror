# Geometry Analysis
## Main Script `main.ipynb`
### Imports & Initialization
Loads estimation and visualization modules; sets core simulation size

#### Variables
+ **Constant case variables**: `big_r` (fixed radius of spheres), `n_sim` (number of simulations), `n_sim_set` - set of simulations repetitions (several `n_sim`), `n_rep`
+ **Variable case variable**s: `theta` & `k` (shape and scale gamma parameter), `n_fol` (number of follicules), `k_set` & `theta_set` (set of gamma parameters)

#### Dependencies
+ **Constant case**: `pdf_const`, `big_r_est`, `est_box_const`
+ **Variable case**: `pdf_gamma`, `viz_gamma_shift`, `viz_gamma_grid`, `grid`, `big_r_gamma_est`, `params_shift`

### Constant‑R Estimation
Visualizes the expected distribution (`pdf_const`), runs a single large simulation (`big_r_est`), and benchmarks performance across varying simulation counts (`est_box_const`)

### Gamma‑R Estimation & Visualization
Explores the impact of gamma parameters on R estimation, generating single‑parameter plots and 2D parameter grids

### Parameter Sweep for Shifted Fits
Systematically fits shifted models varying one parameter at a time, reporting fit parameters across the grid

### How to Run
Consecutive execution of python notebook cells

+ No CLI arguments; adjust parameter values and grids directly in the script before execution
+ Ensure `const_case.py`, `gamma_case.py`, and `viz.py` are in the same directory or on `PYTHONPATH`

# Geometries Dependencies
## Constant Case
**Purpose**: Simulate apparent radii from random planar cross‑sections of a sphere and estimate the true sphere radius

### Overview & Usage
#### Inputs
+ `big_r` (float): true sphere radius
+ `n_sim` (int): number of simulated cross‑sections

#### Outputs
+ `sim(big_r, n_sim)`: array of apparent radii
+ `big_r_est(big_r, n_sim, stats)`: estimated radius (and optional console stats)

### Function Breakdown
#### `sim(big_r, n_sim)`
Draws n_sim distances d from sphere center to slicing plane (uniform in $[0, R]$)

Computes apparent radius $r = \sqrt{R^2 – d^2}$ for each slice

**Returns**: 1D NumPy array of simulated radii

#### `big_r_est(big_r, n_sim, stats = False)`
Uses the mean apparent radius ⟨r⟩ to estimate R via the relation $R ≈ ⟨r⟩·\frac{4}{π}$

If `stats = True`, prints the estimate and absolute error

**Returns**: estimated sphere radius

### How to Use
```py
from radius_estimation import big_r_est
R_true = 5
n_samples = 1000
R_hat = big_r_est(R_true, n_samples, stats = True)
```
+ Adjust `big_r` and `n_sim` to match your sampling requirements

## Variable Case
**Purpose**: Simulate apparent radii from spherical cross‑sections when true radii follow a Gamma distribution, and fit Gamma models to both raw small‐r and back‐estimated large‐R

### Overview & Usage
#### Inputs
+ `theta` (float): shape parameter for true R ~ Gamma(θ, k)
+ `k` (float): scale parameter for true R
+ `n_fol` (int): number of simulated follicles (slices)

#### Outputs
+ `sim_gamma()`: array of small‑r values
+ `gamma_shift()`: fitted Gamma parameters (`a`, `scale`) to small‑r plus raw r_set
+ `big_r_gamma_est()`: prints and returns fitted (`θ_est`, `k_est`) for back‑estimated R

### Function Breakdown
#### `sim_gamma(theta, k, n_fol)`
Draws true radii $R_i ∼ Gamma(\theta, k)$, then for each simulates slice distance $d ∼ Uniform[0, R_i]$, computing apparent radius $r_i = \sqrt{R_i^2 – d^2}$

#### `gamma_shift(theta, k, n_fol)`
Fits a zero‑loc Gamma to the small‑r data, returning estimated shape `a`, `scale`, and the raw `r_set`

#### `big_r_gamma_est(theta, k, n_fol)`
Back‑estimates large radii via $R_{est} = r\bullet \frac{4}{π}$, fits a Gamma to those estimates, and prints estimated vs. true parameters

### How to Use
```py
from gamma_radius_estimation import gamma_shift, big_r_gamma_est
# Example parameters
θ, k, n = 2.0, 3.0, 1000

# Obtain fitted small‑r parameters
a_est, scale_est, r_samples = gamma_shift(θ, k, n)

# Back‑estimate large‑R Gamma parameters and print stats
big_r_gamma_est(θ, k, n)
```

## Visualization Module
**Purpose**: Compare empirical and analytical distributions of apparent sphere‐slice radii, and evaluate Gamma‐based estimators through various visualization routines

### Overview & Usage
#### Inputs
##### Custom modules
+ `const_case`: `sim()`, `big_r_est()`
+ `gamma_case`: `sim_gamma()`, `gamma_shift()`
+ **Modules**: `analytical_pdf()`, `marginal_r()`
+ Simulation parameters hard‑coded at top of the script (e.g., `big_r`, `n_sim_set`, `theta_set`, `k_set`)

#### Outputs
+ PDFs and boxplots saved under `../images/` as: `pdf_const.pdf`, `est_box_const.pdf`, `pdf_gamma.pdf`, `viz_gamma_shift.pdf`, `viz_gamma_grid.pdf`, `params_shift.pdf`

### Pipeline Steps & Key Functions
#### Styling & Theme Setup
Defines a custom color cycle and removes top/right spines for cleaner plots

##### `pdf_const(big_r, n_sim)`
**Simulation**: `r_set = sim(big_r, n_sim)`

**Histogram & Analytical PDF**: plots histogram of `r_set` and overlays `analytical_pdf(r_grid, big_r)`

**Beta Fit**: fits Beta to `r_set`, computes AICs for analytical vs. Beta models, plots fitted Beta PDF

**Save**: writes `../images/pdf_const.pdf`

##### `est_box_const(big_r, n_sim_set, n_rep)`
**Estimation**: runs `big_r_est(big_r, n_sim)` for each `n_sim` in `n_sim_set`, repeated n_rep times

**Boxplots**: produces side‑by‑side boxplots of R estimates for different `n_sim`, marking the true `big_r`

**Save**: writes `../images/est_box_const.pdf`

##### `pdf_gamma(theta, k, n_fol)`
**Simulation**: `r_set = sim_gamma(theta, k, n_fol)`

**Analytical PDF**: computes `marginal_r(r, theta, k)` over a grid

**Plot**: overlays histogram and analytical curve

**Save**: writes `../images/pdf_gamma.pdf`

##### `viz_gamma_shift(theta, k, n_fol, stats = False)`
**Shift Estimation**: uses gamma_shift() to fit a Gamma to r_set

**Plot PDFs**: compares true R‐Gamma vs. fitted r‐Gamma; optionally prints parameter differences

**Save**: writes `../images/viz_gamma_shift.pdf`

##### `viz_gamma_grid(theta_set, k_set, n_fol, mode)`
**Grid Visualization**: for each $(θ, k)$ pair in the 3×3 grid, calls either `pdf_gamma()` or `viz_gamma_shift()` depending on mode

**Save**: writes `../images/viz_gamma_grid.pdf`

##### `params_shift(params_set, mode)`
**Repeated Fitting**: for each parameter value and 30 replicates, runs `gamma_shift()` to gather estimated (`θ`, `k`)

**Aggregation**: averages estimates across replicates

**Plot**: scatter of true vs. estimated values with a y=x reference line

**Save**: writes `../images/params_shift.pdf`

### How to Run
```py
python vis.py # better through main.ipynb
```

+ No command‐line arguments; adjust parameter lists (`n_sim_se`t, `theta_set`, `k_set`, etc.) directly in the script
+ Ensure all imported modules (`const_case.py`, `gamma_case.py`, `modules.py`) are available and dependencies (`seaborn`, `scipy`, `matplotlib`) are installed


