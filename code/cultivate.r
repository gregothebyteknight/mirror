
setwd(this.path::here())

source("nodule.r")
source("vessel.r")

library(ggplot2)

# INITIALIZE VARIABLES
int_nod <- 2000 # intensity of the Poisson process for nodules
int_vess <- 0.005 # intensity of the Poisson process for vessels

num_sim <- 30 # number of simulations
type <- "nodule" # type of simulation ("nodule" or "vessel")
cell_counts <- numeric(num_sim)


# NODULE SIMULATION
for (i in 1:num_sim) {
  dir.create(sprintf("../data/sim_%s", i), showWarnings = FALSE,
             recursive = TRUE)
  if (type == "nodule") { # Simulate nodules
    dot_loc <- nodule(int_nod, viz = FALSE, save = FALSE)
  } else if (type == "vessel") { # Simulate vessels
    dot_loc <- vessel(int_vess, viz = FALSE, save = FALSE)
  }
  cell_counts[i] <- nrow(dot_loc)
  write.csv(dot_loc, sprintf("../data/sim_%s/selected_cell_coordinates.csv", i),
            row.names = FALSE)
}

# VISUALIZATIONs
# Create a data frame for plotting
cell_data <- data.frame(sim_num = 1:num_sim, cell_count = cell_counts)

# Generate the box plot
ggplot(cell_data, aes(x = type, y = cell_count)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(title = "Box Plot of Cell Counts", x = "Simulation Type",
       y = "Number of Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.5))
# Save the plot
ggsave(sprintf("../images/%s_cell_counts_boxplot.png", type),
       width = 8, height = 6)