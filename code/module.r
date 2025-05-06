
# LOADING LIBRARIES
library(spatstat) # for 3D Poisson process
library(Directional) # directional distributions (rkent)

cross <- function(a, b) {
  "compute the cross product of two vectors
  a, b: numeric vectors of length 3"
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

rotation <- function(vec) {
  "compute the rotation matrix for a vector
  with Rodrigues' formula
  vec: numeric vector of length 3"
  axis <- cross(vec, c(1, 0, 0))
  axis_mag <- sqrt(sum(axis^2))
  if (axis_mag < 1e-9) {
    if (vec[1] > 0) return(diag(3))
    else return(diag(c(1, 1, 1)))
  }
  axis <- axis / axis_mag  # Normalize axis
  theta <- acos(sum(vec * c(1, 0, 0)))
  anti_mat <- matrix(c(0, -axis[3], axis[2], axis[3], 0, -axis[1],
                       -axis[2], axis[1], 0), nrow = 3, byrow = TRUE)
  diag(3) + sin(theta) * anti_mat +
    (1 - cos(theta)) * (anti_mat %*% anti_mat)
} # I + sin(theta) * K + (1 - cos(theta)) * K^2 # nolint

pois_in_space <- function(r, d, br_norm, lambda) {
  "generate a 3D Poisson process in a parallelepiped
  r: numeric, radius of the vessel
  d: numeric, thickness of the vessel wall
  br_norm: numeric, length of the direction vector
  lambda: numeric, intensity of the Poisson process"
  box <- spatstat.geom::box3(xrange = c(0, br_norm),
                             yrange = c(-r - d, r + d),
                             zrange = c(-r - d, r + d))
  spatstat.random::rpoispp3(lambda, domain = box)
}

kent_point <- function(n, kappa, vec, beta, last_dot, step_size) {
  "generate a new point using the Kent distribution
  n: integer, number of points to generate
  kappa: numeric, concentration parameter,
    higher values result in more concentrated points
  vec: numeric vector of length 3, direction vector
  beta: numeric, ovalness parameter, 
    higher values result in more elongated points
  last_dot: numeric vector of length 3, last point
  step_size: numeric, distance to extend the branch"
  vec <- Directional::rkent(n, kappa, vec / sqrt(sum(vec^2)), beta)[1, ]
  last_dot + vec * step_size
}