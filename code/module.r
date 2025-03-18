
# LOADING LIBRARIES
library(pracma) # for cross function

deg2rad <- function(deg) (deg * pi) / (180)

cross <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

rotation <- function(vec) { # Rodrigues' rotation formula
  axis <- pracma::cross(vec, c(1, 0, 0))
  axis_mag <- sqrt(sum(axis^2))
  if (axis_mag < 1e-9) {
    if (vec[1] > 0) return(diag(3))
    else return(diag(c(-1, -1, 1)))
  }
  axis <- axis / axis_mag  # Normalize axis
  theta <- acos(sum(vec * c(1, 0, 0)))
  anti_mat <- matrix(c(0, -axis[3], axis[2], axis[3], 0, -axis[1],
                       -axis[2], axis[1], 0), nrow = 3, byrow = TRUE)
  diag(3) + sin(theta) * anti_mat +
    (1 - cos(theta)) * (anti_mat %*% anti_mat)
} # I + sin(theta) * K + (1 - cos(theta)) * K^2 # nolint

# pois_space <- function(w) {
#   box <- spatstat.geom::box3(
#     xrange = c(0, w),
#     yrange = c(0, 2 * r),
#     zrange = c(0, 2 * r)
#   )
#     prizm <- rpoispp3(lambda, domain = box)
#     prizm <- cbind(prizm$data$x, prizm$data$y - r, prizm$data$z - r)  # Center
# }