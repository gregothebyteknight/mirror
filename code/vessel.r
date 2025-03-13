
# LOADING LIBRARIES
library(rgl)         # interactive 3D plotting
library(fifer)       # for scring.to.colors
library(Directional) # directional distributions (rkent)
library(htmlwidgets)

setwd(this.path::here())
options(rgl.printRglwidget = TRUE)

# VARIABLES DECLARATION
t_total <- 10
delta_t <- 0.01
step_size <- 1  # step size
prob <- 0.005
num_steps <- as.integer(t_total / delta_t)

vec_init <- rkent(2, 0, c(1, 0, 0), 1)[1, ]  # initial direction
tree <- list()                        # list to hold branch trajectories
tree[[1]] <- matrix(c(0, 0, 0), nrow = 1)  # start with origin
tree[[1]] <- rbind(tree[[1]], vec_init * step_size)

for (step in 1:num_steps) {
  new_branches <- list()
  current_tree_length <- length(tree)
  # EXTENDING EXISTING BRANCHES
  for (i in 1:current_tree_length) {
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
dot_list_total <- do.call(rbind, tree)
# Create a batch vector that assigns a unique id to each branch’s points
batch <- unlist(lapply(seq_along(tree), function(i) rep(i, nrow(tree[[i]]))))

plot3d(dot_list_total, col = string.to.colors(batch))
# Save as HTML widget
widget <- rglwidget()
saveWidget(widget, "../images/vessel.html", selfcontained = TRUE)
