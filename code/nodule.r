library(rgl)
library(htmlwidgets)

setwd(this.path::here())
options(rgl.printRglwidget = TRUE)

# INITIALIZE VARIABLES
border <- 1000
win_min <- 0
win_max <- 1000
win_scale <- win_max - win_min
inits <- rpois(n = 1, lambda = 30) # number of clusters
init_loc <- cbind(runif(inits, win_min, win_max),
                  runif(inits, win_min, win_max),
                  runif(inits, win_min, win_max)) # location of inits
dot_loc <- matrix(ncol = 3, nrow = 0) # locations of dots

# GENERATE SPHERES OF DOTS AROUND INITIALES
for (init in 1:inits) {
  r <- rexp(n = 1, rate = 1 / (0.1 * win_scale)) # random radius
  dot_temp <- rpois(n = 1, lambda = (4 / 3) * pi * r^3 * 5000 / win_scale^3)

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
rgl::plot3d(dot_loc, col = "#5454c4e1", size = 2,
            xlab = "X", ylab = "Y", zlab = "Z")

# Save as HTML widget
widget <- rgl::rglwidget()
htmlwidgets::saveWidget(widget, "../images/nodule.html", selfcontained = TRUE)

# Save as CSV
dot_loc <- scale(dot_loc, center = TRUE, scale = FALSE) # centering
colnames(dot_loc) <- c("x", "y", "z") # rename columns
utils::write.csv(dot_loc, "../data/nodule_loc.csv", row.names = FALSE)