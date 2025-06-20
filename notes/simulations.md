# Master Script
**Purpose**: Simulate spatial point patterns of either “nodules” or “vessels” via Poisson processes, save each simulation’s cell coordinates, and summarize counts in a box plot

## Overview & Usage
### Inputs (hard‑coded)
+ `int_nod = 7e4`: intensity (points per unit volume) for nodules
+ `int_vess = 5e-3`: intensity for vessels
+ `num_sim = 30`: number of replicate simulations
+ type = "nodule" or "vessel": choose simulation model

### Outputs
+ Per-simulation CSVs in `../data/sim_<i>/selected_cell_coordinates.csv`
+ Boxplot PNG at `../images/<type>_cell_counts_boxplot.png`

## Pipeline Steps & Key Logic
### Setup & Imports
Sets project root, loads the two simulation modules, and the plotting library

### Parameter Initialization
Defines process intensities, number of replicates, and preallocates a vector to record simulated point counts

### Simulation Loop
For each replicate: creates output folder, runs the appropriate point‐pattern simulation, records the number of points, and saves the coordinate table

### Summarize & Plot
Constructs a simple box plot comparing replicate counts, then saves it as a PNG

## How to Run
```bash
Rscript simulate_and_plot.r
```
+ Ensure nodule.r and vessel.r are present (dependencies)
+ Adjust `type`, `int_nod`/`int_ves`, and `num_sim` at the top to explore other scenarios

# Nodules `nodule.r`
**Purpose**: Generate clustered 3D “nodule” point patterns via a Poisson cluster process, optionally visualize in an interactive HTML widget, and center/save coordinates with cluster labels

## Overview & Usage
### Inputs (function arguments)
+ `int` (numeric): overall point intensity (points per unit³)
+ `ord` (integer): simulation index, used in the returned cell label
+ `viz` (logical): if TRUE, render asnd save a WebGL HTML of the 3D points
+ `save` (logical): if TRUE, write centered coordinates to `../data/nodule_loc.csv`

### Output
+ Invisibly returns a data matrix (when save = FALSE) with columns `x`, `y`, `z`, `cluster`, `cell`
+ When `viz = TRUE`, creates `../images/nodule.html`
+ When `save = TRUE`, writes `../data/nodule_loc.csv`

## Function Breakdown: `nodule(int, ord, viz = FALSE, save = FALSE)`
### Window & Cluster Initialization
Defines a cubic window $[0,500]^3$, simulates Poisson cluster centers uniformly in the volume

### Clustered Point Generation
For each cluster center, draws an Exponential radius, then populates that sphere with a Poisson number of points, accepting only those inside the sphere

### Add Noise & Center
Adds Gaussian jitter (5% of window size) and recenters coordinates to zero-mean

### Visualization (Optional)
Launches an interactive 3D widget via `rgl` and saves it as a self‑contained HTML

### Label & Save
Appends cluster and cell labels, then writes to CSV if requested; otherwise returns the data

## How to Run
```r
source("nodule_simulation.r")
nodule(int = 7e4, ord = 1, viz = TRUE, save = TRUE)
```

+ Ensure `rgl` and `htmlwidgets` are installed for 3D visualization
+ Adjust `int` or `lambda` in `rpois()` for different cluster counts

# Vessels `vessel.r`
**Purpose**: Simulate 3D “vessel” point patterns by growing branching trajectories and populating vessel walls with Poisson‐distributed points; optionally visualize in WebGL and save coordinates

## Overview & Usage
### Inputs (function arguments)
+ `int` (numeric): Poisson intensity for wall‐points per unit volume
+ `ord` (integer): simulation identifier for labeling
+ `viz` (logical): if TRUE, render an interactive 3D widget and save as HTML
+ `save` (logical): if TRUE, write centered coordinates to `../data/vessel_loc.csv`

### Output
+ Invisibly returns a data matrix with columns `x`, `y`, `z`, `cluster`, `cell`
+ If `viz = TRUE`, creates `../images/vessel_combined.html`
+ If `save = TRUE`, writes `../data/vessel_loc.csv`

## Pipeline Steps & Key Logic
### Setup & Imports
Establish project root, load dependence on custom modules (`pois_in_space`, `rotation`, etc.), 3D plotting and Poisson‐process utilities

### Branch Trajectory Growth
Builds a list of branches, each as a matrix of 3D points
```r
tree <- list(matrix(c(0, 0, 0), 1))
for(step in 1:200) {
  for(each branch in tree) {
    extend branch by one step via `kent_point()`
    with probability 0.006, create a new branch off the current tip
  }
}
```

### Wall‐Point Generation
For each linear piece of every branch, simulates a Poisson parallelepiped, retains points near the vessel wall, and positions them in 3D
```r
    for(each branch segment) {
      prizm <- pois_in_space(r = 7.5, d = 1.5, length = segment_length, lambda = int)
      keep points within wall thickness ±0.75
      rotate points to align with segment via `rotation()`
      translate to segment start and add Gaussian jitter (SD = 0.5)
      append to `dot_loc`
    }
```

### Visualization (Optional)
Renders both the vessel skeleton and wall‐points in an interactive widget
```r
if(viz) {
  plot3d(branches, col=string.to.colors(batch))
  points3d(dot_loc)
  saveWidget(rglwidget(), "../images/vessel_combined.html")
}
```

### Centering & Labeling
Recenters coordinates to zero‐mean, adds annotation columns, and optionally saves CSV
```r
    dot_loc <- scale(dot_loc, center=TRUE, scale=FALSE)
    dot_loc <- cbind(dot_loc, cluster=1, cell=sprintf("Vessel cell_%s", ord))
    colnames(dot_loc) <- c("x","y","z","cluster","cell")
    if(save) write.csv(dot_loc, "../data/vessel_loc.csv", row.names=FALSE)
```

## How to Run
```r
source("vessel.r")
vessel(int = 5e-3, ord = 1, viz = TRUE, save = FALSE)
```
+ Ensure `kent_point()`, `pois_in_space()`, and `rotation()` are defined in `module.r`
+ Install required packages (`rgl`, `fifer`, `htmlwidgets`, `spatstat`) for full functionality