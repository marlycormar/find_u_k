
install.packages("primes")
install.packages("rlist")
library(primes)
library(MASS)
library(rlist)

# Initial setup
max_prime <- 200 # max prime to consider
primes <- primes::generate_primes(2, max_prime) # get the first max_prime primes
p_1 <- primes[1] # get the first prime

# First elasticities
u_1 <- c(1)
elasticities <- list(u_1)

# Min element in U_k
min_elast <- function(k){
  return (ceiling(k*(p_1 - 1)/p_1))
}

# Max element in U_k
max_elast <- function(k){
  return (floor(k*p_1/(p_1 - 1)))
}

# Testing
#min_elast(5)
#max_elast(5)


# Traverses the elements l in U_{k-1} to add l+1 to U_k, and 
# k-1+diff to U_k, where diff=k-l
add_known_elasticities_from_previous_sets <- function(k) {
  #k <- 5
  u_k <- list()
  u <- unlist(elasticities[k-1])
  
  for(j in 1:length(u)){
    #j <- 2
    u_k[length(u_k) + 1] <- u[j] + 1
    
    if(u[j] <= k){
      u_k[length(u_k) + 1] <- k-1 + k-u[j]
    }
  }

  return(unique(u_k))
}

# Testing
#u_5 <- add_known_elasticities_from_previous_sets(5)

# Since we are only considering when n is an integer
# we can find the max. atom that can be in a factorization 
# of n, namely (p_i-1)/p_i where i is the smallest index 
# such that p_i-1 >= n. 
max_atom_index_to_consider <- function(n){
  #n <- 4
  for(i in 1:n){
    #i <- 4
    if((primes[i] - 1) >= n){
      return(i)
    }
  }
  print("It didnt find the max_atom_to_consider.")
  return(0)
}


# The list atoms contains only the numerators of atoms
length_of_n_as_the_sum_of_elements_in_list_rec <- function(n, max_atom_index_to_consider, sum, total_summands){
  if(n == sum) return (c(total_summands))
  if(sum > n) return (c(NA))
  result <- list()
  for(i in 1:max_atom_index_to_consider){
    temp <- length_of_n_as_the_sum_of_elements_in_list_rec(n, max_atom_index_to_consider, sum + primes[i] - 1, total_summands + primes[i])
    result[[length(result) + 1]] <- temp
  }
  return (result)
}


length_of_n_as_the_sum_of_elements_in_list <- function(n, max_atom_index_to_consider){
  if(max_atom_index_to_consider == 0) return (c(NA))
  if(n < (primes[1] - 1)) return (c(NA))
  
  res <- unlist(length_of_n_as_the_sum_of_elements_in_list_rec(n, max_atom_index_to_consider, 0, 0))
  res <- unique(res[!is.na(res)])
  
  return (res)
}

# Testing
#res <- length_of_n_as_the_sum_of_elements_in_list_rec(4, c(1, 2, 4), 0, 0)
#res <- length_of_n_as_the_sum_of_elements_in_list(4, c(1, 2, 4))

add_new_elasticities <- function(k){
  #k <- 5
  #n <- 4
  u_k <- list()
  #possible_elast <- setdiff(min_elast(k):max_elast(k), add_known_elasticities_from_previous_sets(k))
  max_integer <- k - 1
  min_integer <- ceiling(k/2)
  # Each atom must be multiplied by their denominator (check this)
  for(n in min_integer:max_integer){
    max_atom_index_to_consider <- max_atom_index_to_consider(n)
    u_k <- c(u_k, setdiff(length_of_n_as_the_sum_of_elements_in_list(n, max_atom_index_to_consider), u_k))
  }
  return (u_k)
}

find_u_k <- function(k){
  u_k <- unique(unlist(c(add_known_elasticities_from_previous_sets(k), add_new_elasticities(k))))
  return (sort(u_k))
}

draw_plot <- function(){
  dim <- 2*length(elasticities)
  plot(1:dim, 1:dim, type = "n")  # setting up coord. system
  for(i in 1:length(elasticities))
  {
    for(j in 1:length(elasticities[[i]]))
    {
      points(i, (elasticities[[i]])[j], col = "red", pch=19)
    }
  }
}

min_k <- 24 # from which elasticities we will start
max_k <- 25
for(i in min_k:max_k){
  elasticities[[i]] <- find_u_k(i)
}
names(elasticities) <- 1:length(elasticities)
elasticities <- `20_first_elasticities`
draw_plot()

#saveRDS(elasticities, file = "20_first_elasticities.rds")


# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

length(elasticities[[i]])

dim <- 2*length(elasticities)
plot(1:dim, 1:dim, type = "n")
i <- 3
j <- 1
points(i, (elasticities[[i]])[j], col = "red", pch=19)
i <- 3
j <- 2
points(i, (elasticities[[i]])[j], col = "red", pch=19)


# Crazy testing
z <- 1:5
c <- c(1, 5, 10)
setdiff(z, c)
setdiff(c, z)
str(z)
c <- z[!(z %in% c(2, 3))]
c

temp <- list()
other <- list(2, 3, 4)
temp[[1]] <- other
temp[[2]] <- c(other, 100)
unlist(temp)


sort(c(4,3))

# Simple Scatterplot
# Simple Scatterplot
attach(mtcars)
plot(wt, mpg, main="Scatterplot Example", 
  	xlab="Car Weight ", ylab="Miles Per Gallon ", pch=19)

names(elasticities) <- 1:length(elasticities)

