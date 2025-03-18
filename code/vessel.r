
# LOADING LIBRARIES
source("./module.r")

library(rgl)         # interactive 3D plotting
library(fifer)       # for scring.to.colors
library(spatstat)
library(data.table)
library(Directional) # directional distributions (rkent)
library(htmlwidgets)

setwd(this.path::here())
options(rgl.printRglwidget = TRUE)

# VARIABLES DECLARATION
t_total <- 10
prob <- 0.005
step_size <- 0.5  # step size
delta_t <- 0.025
num_steps <- as.integer(t_total / delta_t)

r <- 3 # radius of vessel
d <- 1 # width of vessel's walls
lambda <- 0.02 # for Poisson process
jitter_sd <- 0.1  # standard deviation for jitter

vec_init <- rkent(2, 0, c(1, 0, 0), 1)[1, ]  # initial direction
tree <- list() # list to hold branch trajectories
tree[[1]] <- matrix(c(0, 0, 0), nrow = 1)  # start with origin
tree[[1]] <- rbind(tree[[1]], vec_init * step_size)

# INTERATION THROUGH THE TREE
for (step in 1:num_steps) {
  new_branches <- list()
  tree_len <- length(tree)
  # EXTENDING EXISTING BRANCHES
  for (i in 1:tree_len) {
    branch <- tree[[i]]
    last_dot <- tree[[i]][nrow(branch), ]
    # Compute the current direction from the last two points (if available)
    if (nrow(tree[[i]]) > 1) {
      vec <- last_dot - tree[[i]][nrow(tree[[i]]) - 1, ]
      vec <- vec / sqrt(sum(vec^2))
    } else {
      vec <- vec_init
    }
    # Extend the branch using a new direction based on the current direction
    new_dir <- rkent(50, 1000, vec, 0)[1, ]
    new_point <- last_dot + new_dir * step_size
    tree[[i]] <- rbind(tree[[i]], new_point)

    # BRANCHING CONDITION
    # Branching event: create a new branch starting from the current tip
    if (runif(1) < prob) {
      new_branch_dir <- rkent(2, 0, vec, 1)[1, ]
      branch_new_point <- last_dot + new_branch_dir * step_size
      branch <- rbind(branch, branch_new_point)
      new_branches[[length(new_branches) + 1]] <- branch
    }
  }
  tree <- c(tree, new_branches)
}

# VISUALIZATION
# Combine all branch points for plotting
init_dot_loc <- do.call(rbind, tree) # set of dots forming axis of vessels
# Create a batch vector that assigns a unique id to each branch’s points
batch <- unlist(lapply(seq_along(tree), function(i) rep(i, nrow(tree[[i]]))))

rgl::plot3d(init_dot_loc, col = string.to.colors(batch),
            xlab = "X", ylab = "Y", zlab = "Z")
# Save as HTML widget
widget <- rgl::rglwidget()
htmlwidgets::saveWidget(widget, "../images/vessel.html", selfcontained = TRUE)

# VESSEL SIMULATION
dot_loc <- matrix(ncol = 3, nrow = 0) # final matrix of dot location

for (branch in tree) {
  for (i in seq_len(nrow(branch) - 1)) {
    # Compute branch direction vector
    curr_pnt <- branch[i, ]
    next_pnt <- branch[i + 1, ]
    br_vec <- next_pnt - curr_pnt
    br_norm <- sqrt(sum(br_vec^2))
    br_vec <- br_vec / br_norm

    # modify sample poisson parallelepiped
    # Generate 3D Poisson parallelepiped
    box <- spatstat.geom::box3(
      xrange = c(0, br_norm),
      yrange = c(-r - d, r + d),
      zrange = c(-r - d, r + d)
    )

    prizm <- spatstat.random::rpoispp3(lambda, domain = box)
    if ((npoints(prizm) > 0)) {  # Check if any points are generated
      prizm <- cbind(
        as.numeric(prizm$data$x),
        as.numeric(prizm$data$y),
        as.numeric(prizm$data$z)
      )

      # Modify to tube
      cond <- sqrt(prizm[, 2]^2 + prizm[, 3]^2)
      prizm <- prizm[abs(cond - r) <= d / 2, , drop = FALSE]

      # Compute and apply rotation
      prizm_rot <- t(rotation(br_vec) %*% t(prizm))

      # Translation
      prizm_rot[, 1] <- prizm_rot[, 1] + curr_pnt[1]
      prizm_rot[, 2] <- prizm_rot[, 2] + curr_pnt[2]
      prizm_rot[, 3] <- prizm_rot[, 3] + curr_pnt[3]
      prizm_rot <- prizm_rot +
        matrix(rnorm(length(prizm_rot), 0, jitter_sd), ncol = 3)

      dot_loc <- rbind(dot_loc, prizm_rot)
    }
  }
}

rgl::plot3d(dot_loc, xlab = "X", ylab = "Y", zlab = "Z", size = 2)
widget <- rgl::rglwidget()
htmlwidgets::saveWidget(widget, "../images/vessel_dots.html",
                        selfcontained = TRUE)