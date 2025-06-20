
# SETUP THE ENVIRONMENT
setwd(this.path::here())
suppressPackageStartupMessages({
  library(rgl)
  library(htmlwidgets)
})
options(rgl.printRglwidget = TRUE)

nodule <- function(int, ord, viz = FALSE, save = FALSE) {
  # VARIABLES DECLARATION
  win_min <- 0
  win_max <- 500
  win_scale <- win_max - win_min
  inits <- rpois(n = 1, lambda = 15) # number of clusters
  init_loc <- cbind(runif(inits, win_min, win_max),
                    runif(inits, win_min, win_max),
                    runif(inits, win_min, win_max)) # location of inits
  dot_loc <- matrix(ncol = 3, nrow = 0) # locations of dots

  # GENERATE SPHERES OF DOTS AROUND INITIALES
  for (init in 1:inits) {
    r <- rexp(n = 1, rate = 1 / (0.05 * win_scale)) # random radius
    dot_temp <- rpois(n = 1, lambda = (4 / 3) * pi * r^3 * int / win_scale^3)

    i <- 0
    while (i < dot_temp) {
      loc_temp <- runif(n = 3, min = init_loc[init, ] - r,
                        max = init_loc[init, ] + r)
      dot_r_temp <- sqrt(sum((loc_temp - init_loc[init, ])^2))
      if (dot_r_temp <= r) { # euclid dist <= radius accepted
        dot_loc <- rbind(dot_loc, loc_temp)
        i <- i + 1
      }
    }
  }

  dot_loc <- dot_loc + matrix(data = rnorm(n = nrow(dot_loc) * 3,
                                           mean = 0, sd = 0.05 * win_scale),
                              nrow = nrow(dot_loc), ncol = 3) # add noise
  # VISUALIZE
  if (viz) {
    rgl::plot3d(dot_loc, col = "#5454c4e1", size = 2,
                xlab = "X", ylab = "Y", zlab = "Z")

    # Save as HTML widget
    widget <- rgl::rglwidget()
    htmlwidgets::saveWidget(widget, "../images/nodule.html",
                            selfcontained = TRUE)
  }

  # Final edits
  dot_loc <- scale(dot_loc, center = TRUE, scale = FALSE) # centering
  dot_loc <- cbind(dot_loc, rep("1", nrow(dot_loc))) # add clusters
  dot_loc <- cbind(dot_loc, rep(sprintf("Nodule cell_%s", ord), nrow(dot_loc)))
  colnames(dot_loc) <- c("x", "y", "z", "cluster", "cell") # rename columns
  print(sprintf("Number of cells in the simulation: %s", nrow(dot_loc)))
  if (save) {
    utils::write.csv(dot_loc, "../data/nodule_loc.csv", row.names = FALSE)
  } else {
    return(invisible(dot_loc))
  }
}


nodule(7e4, 1, viz = TRUE, save = TRUE)
