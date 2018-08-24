install.packages("primes")
install.packages("rlist")
library(primes)
library(MASS)
library(rlist)


# Variables to modify
smallest_prime_to_consider <- 2
max_k <- 20 # Max k for which to find U_k

# Initial setup
max_prime <- 200 # max prime to consider
primes <- primes::generate_primes(smallest_prime_to_consider, max_prime) # get the first max_prime primes
p_1 <- primes[1] # get the first prime

# Traverses the elements l in U_{k-1} to add l+1 to U_k, and 
# k-1+diff to U_k, where diff=k-l
# u_k_minus_1 = U_{k-1}
add_known_lengths_from_u_k_minus_1 <- function(k, u_k_minus_1) {
  #k <- 2
  u_k <- list()
  for(j in 1:length(u_k_minus_1)){
    #j <- 1
    u_k[length(u_k) + 1] <- u_k_minus_1[j] + 1
    
    if(u_k_minus_1[j] <= k && u_k_minus_1[j] != (k - 1)){
      u_k[length(u_k) + 1] <- k-1 + k-u_k_minus_1[j]
    }
  }

  return(unique(u_k))
}

# Testing
#u_4 <- add_known_lengths_from_u_k_minus_1(4, c(3, 4))
#u_3 <- add_known_lengths_from_u_k_minus_1(3, c(2))

# Since we are only considering the case when n is an integer
# we can find the max. atom that can be in a factorization 
# of n, namely (p_i-1)/p_i where i is the smallest index 
# such that p_i-1 >= n. 
# Returns the index of such atom.
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

# n: an integer that is assumed to have a factorization of length k
#    and the one we are finding all possible lengths of factorizations of.
# max_atom_index_to_consider: the index of the max. possible atom that can be in a factorization of n.
# sum: current factorization sum.
# fact_length: current factorization length.
add_new_lengths_for_n_rec <- function(n, max_atom_index_to_consider, sum, fact_length){
  if(n == sum) return (c(fact_length))
  if(sum > n) return (c(NA))
  result <- list()
  for(i in 1:max_atom_index_to_consider){
    temp <- add_new_lengths_for_n_rec(n, max_atom_index_to_consider, sum + primes[i] - 1, fact_length + primes[i])
    result[[length(result) + 1]] <- temp
  }
  return (result)
}

# Adds the lengths not found with method add_known_lengths_from_u_k_minus_1.
# By our observations, if $x$ has a factorization of orden $k$ not obtained from $U_{k-1}$,
# then $x$ must be an integer in certain interval. So it sufficies to iterate through such
# interval finding all the possible lenghts of elements in that interval.
add_new_lengths_to_u_k <- function(k){
  #k <- 5
  #n <- 4
  u_k <- list()
  max_integer <- k - 1 # See observations.
  min_integer <- ceiling(k*(p_1 - 1)/p_1) # See observations.
  # Each atom must be multiplied by a multiple of its denominator.
  for(n in min_integer:max_integer){
    max_atom_index_to_consider <- max_atom_index_to_consider(n)
    if(max_atom_index_to_consider == 0 || (n < (primes[1] - 1))){# Do some basic checks
      print("Either max_atom_index_to_consider = 0 or n < (primes[1] - 1)")
    }
    else{
      res <- unlist(add_new_lengths_for_n_rec(n, max_atom_index_to_consider, 0, 0))
      res <- unique(res[!is.na(res)]) # Remove NA entries.
      if(k %in% res)
        u_k <- c(u_k, setdiff(res, u_k)) # Add only the new values of integers that contain a factorization of length k.
    }
  }
  return (u_k)
}

# Finds u_k given u_{k-1}
find_u_k <- function(k, u_k_minus_1){
  u_k <- unique(unlist(c(add_known_lengths_from_u_k_minus_1(k, u_k_minus_1), add_new_lengths_to_u_k(k))))
  return (sort(u_k))
}

# Testing
# <- find_u_k(3, c(2))


# Finds u_k for k between 1 and max_k.
# Returns a named list containing {u_1, u_2, }
find_set_of_union_of_sets_of_lengths <- function(max_k){
  #max_k = 3
  set_of_unions <- list()
  set_of_unions[1] <- list(1)
  for(i in 2:max_k){
    set_of_unions[[i]] <- find_u_k(i, set_of_unions[[i-1]])
  }
  names(set_of_unions) <- 1:max_k
  return(set_of_unions)
}

find_set_of_union_of_sets_of_lengths_from_min_k_to_max_k <- function(min_k, max_k, set_of_union_of_sets_of_lengths_for_k_less_than_min_k){
  #max_k = 3
  for(i in min_k:max_k){
    set_of_union_of_sets_of_lengths_for_k_less_than_min_k[[i]] <- find_u_k(i, set_of_union_of_sets_of_lengths_for_k_less_than_min_k[[i-1]])
  }
  names(set_of_union_of_sets_of_lengths_for_k_less_than_min_k) <- 1:max_k
  return(set_of_union_of_sets_of_lengths_for_k_less_than_min_k)
}

draw_plot <- function(set_of_union_of_sets_of_lengths){
  height <- 2*length(set_of_union_of_sets_of_lengths)
  width <- length(set_of_union_of_sets_of_lengths)
  xValues <- c()
  yValues <- c()
  dotsize <- 1

  for(i in 1:length(set_of_union_of_sets_of_lengths))
  {
    for(j in 1:length(set_of_union_of_sets_of_lengths[[i]]))
    {
      xValues <- c(xValues, i)
      yValues <- c(yValues, (set_of_union_of_sets_of_lengths[[i]])[j])
    }
  }
  
  dev.new(height+10,width+20)
  plot(xValues, yValues, cex=dotsize, xlim=c(0,width+1), ylim=c(0, height+1), bty="n", 
       pch=18, col=2,  xlab="k", ylab="Uk", axes = FALSE)
  #plot(xValues, yValues, cex=dotsize, xlim=c(0,width+1), ylim=c(0, height+1), 
  #     bty="n", pch=1,  col=cm.colors(2), xlab="k", ylab="Uk", axes = FALSE)
  axis(side = 1, at = 1:width)
  #box()
  axis(side = 2, at = 1:height)
}

set_of_union_of_sets_of_lengths <- c()
set_of_union_of_sets_of_lengths <- find_set_of_union_of_sets_of_lengths(max_k)
draw_plot(set_of_union_of_sets_of_lengths)

set_of_union_of_sets_of_lengths_2 <- find_set_of_union_of_sets_of_lengths_from_min_k_to_max_k(21, 22, set_of_union_of_sets_of_lengths)
draw_plot(set_of_union_of_sets_of_lengths_2)

#saveRDS(set_of_union_of_sets_of_lengths, file = "from_u1_to_u20.rds")
#saveRDS(set_of_union_of_sets_of_lengths_2, file = "from_u21_to_u22.rds")