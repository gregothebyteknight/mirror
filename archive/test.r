

library(rgl)
library(pracma)
library(spatstat)
options(rgl.printRglwidget = TRUE)
lambda <- 1

r <- 5
dot_loc <- matrix(ncol = 3, nrow = 0) # final matrix of dot location

box <- spatstat.geom::box3(xrange = c(0, 2 * r), yrange = c(0, 2 * r),
                           zrange = c(0, 10))

pts <- rpoispp3(lambda, domain = box)

plot3d(pts$data, xlab = "X", ylab = "Y", zlab = "Z", size = 2)
widget <- rglwidget()

pts_coords <- cbind(pts$data$x, pts$data$y, pts$data$z)
head(pts_coords)