library(tidyverse)
library(vegan)

data <- read.csv('Data_JEcol.csv', sep = ';') # loading the data stored here https://doi.org/10.5061/dryad.8gtht76r9 

##### data manipulation #####
data <- 
  data |>
  mutate(across(starts_with('X'), ~ log10(1000 * .x + 1))) # Here with log transformed data

data_long <- # not necessary, I was using this in other analyses 
  data |>
  pivot_longer(cols = starts_with('X'), names_to = 'Year', values_to = 'Biomass') |>
  mutate(Year = str_remove(Year, "^X") |> as.numeric())

com_data <- data_long |> # this is the relevant input for the analyses
  group_by(Plot, Year, Species) |>
  summarise(Biomass = sum(Biomass), .groups = "drop") |>
  pivot_wider(names_from = Species, values_from = Biomass, values_fill = 0)

head(com_data)

##### Two term locl quadrat variance function for later use #####
calc_2tlqv <- function(x) {
  n <- length(x)
  
  diffs <- (x[1:(n - 1)] - x[2:n])^2
  
  sum(diffs) / (2 * (n - 1))
}


##### Euclidean CV #####

plot_ids <- c(1, 3, 6, 8, 9, 11, 14, 16, 18, 20, 21, 23) # List of plots to analyze

results_euclid <- map_dfr(plot_ids, 
                          function(i) { # going across individual plots
  plot_data <- com_data |> # taking the community matrix (log transformed data here)
    filter(Plot == i) |> # selecting one plot
    select(-c(Year, Plot)) |> # keeping just the community mtrix
    as.matrix()
  
  mean_vec <- colMeans(plot_data) # Centroids of species in time
  
  data_merged <- rbind(mean_vec, plot_data) # adding centroid to community matrix
  
  dist_matrix <- as.matrix(vegdist(data_merged, method = "euclidean")) # calculating distance matrix
  dist_euclid <- dist_matrix[-1, 1] # saving just the distance of each point from centroid
  
  var_euclid <- sum(dist_euclid^2) / (nrow(plot_data) - 1) # euclidean variance
  cv_euclid <- sqrt(var_euclid) / sqrt(sum(mean_vec^2)) # Dividing by sample norm
  
  tibble(
    Plot = i,
    comp_var_euclid = var_euclid,
    euclid_CV = cv_euclid,
    sample_norm = sqrt(sum(mean_vec^2)),
    dist_matrix = list(dist_matrix) # distance matrix needed for the next step - i.e. "comtrendr" by 2tlqv
  )
}
)

head(results_euclid)

##### Two term variance for Euclid #####
# getting superdiagonal - so we get between-year euclid distances
superdiagonals <- map(results_euclid$dist_matrix, function(mat) {
  mat <- mat[-1, -1] # getting rid of the centroid again
  diag(mat[, -1]) # deleting one column so superdiagonal becomes diagonal
})

euclid_2tlqv <- data.frame(euclid_2tlqv = unlist(lapply(superdiagonals, calc_2tlqv)), # function for two term variance from above
                           Plot = c(1, 3, 6, 8, 9, 11, 14, 16, 18, 20, 21, 23),
                           sample_norm = results_euclid$sample_norm)

euclid_2tlqv <- 
  euclid_2tlqv %>%
  mutate(euclid_2tlqv_CV = sqrt(euclid_2tlqv) / sample_norm) # CV of 2tlqv

euclid_2tlqv


##### Chord distance #####
# almost the same as for euclidean, I comment just on the differences
results_chord <- map_dfr(plot_ids, function(i) {
  plot_data <- com_data %>%
    filter(Plot == i) %>%
    select(-c(Year, Plot)) %>%
    as.matrix()
  
  mean_vec <- colMeans(plot_data)
  
  data_merged <- rbind(mean_vec, plot_data)
  
  dist_matrix <- as.matrix(vegdist(data_merged, method = "chord")) # calculating chord distance
  dist_chord <- dist_matrix[-1, 1]
  
  var_chord <- sum(dist_chord^2) / (nrow(plot_data) - 1)
  cv_chord <- sqrt(var_chord)  # chord distances are already normalized so no need to divide by sample norm
  
  tibble(
    Plot = i,
    comp_var_chord = var_chord,
    chord_CV = cv_chord,
    dist_matrix = list(dist_matrix)
  )
})

head(results_chord)


##### Two term variance for chord #####
superdiagonals <- map(results_chord$dist_matrix, function(mat) {
  mat <- mat[-1, -1]
  diag(mat[, -1]) 
})


chord_2tlqv <- data.frame(chord_2tlqv = unlist(lapply(superdiagonals, calc_2tlqv)),
                          Plot = c(1, 3, 6, 8, 9, 11, 14, 16, 18, 20, 21, 23))

head(chord_2tlqv)

