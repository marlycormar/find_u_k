# Add the R script containing the functions used here 
source("functions.R")

# Variables to modify
smallest_prime_to_consider <- 2
max_k <- 20 # Max k for which to find U_k
max_prime <- 200 # max prime to consider

# Find U_k for each k in [1:max_k] and make a plot
set_of_union_of_sets_of_lengths <- find_set_of_union_of_sets_of_lengths(max_k)
print(set_of_union_of_sets_of_lengths)
draw_plot(set_of_union_of_sets_of_lengths)
