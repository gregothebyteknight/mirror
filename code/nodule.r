library(rgl)
library(htmlwidgets)
setwd(this.path::here())

options(rgl.printRglwidget = TRUE)

# INITIALIZE VARIABLES
inits <- rpois(n = 1, lambda = 30) # number of clusters
init_loc <- cbind(runif(inits, 0, 1), runif(inits, 0, 1),
                  runif(inits, 0, 1)) # location of inits
dot_loc <- matrix(ncol = 3, nrow = 0) # locations of dots

# GENERATE SPHERES OF DOTS AROUND INITIALES
for (init in 1:inits) {
  r <- rexp(n = 1, rate = 1 / 0.1) # random radius
  dot_temp <- rpois(n = 1, lambda = r^2 * pi * 5000)

  i <- 0
  while (i < dot_temp) {
    loc_temp <- runif(n = 3, min = 0, max = 1)
    dot_r_temp <- sqrt(sum((loc_temp - init_loc[init, ])^2))
    if (dot_r_temp <= r) { # euclid dist <= radius accepted
      dot_loc <- rbind(dot_loc, loc_temp)
      i <- i + 1
    }
  }
}

dot_loc <- dot_loc + matrix(data = rnorm(n = nrow(dot_loc) * 3,
                                         mean = 0, sd = 0.05),
                            nrow = nrow(dot_loc), ncol = 3)
plot3d(dot_loc, col = "#5454c4e1", size = 2,
       xlab = "X", ylab = "Y", zlab = "Z")

# Save as HTML widget
widget <- rglwidget()
saveWidget(widget, "../images/nodule.html", selfcontained = TRUE)