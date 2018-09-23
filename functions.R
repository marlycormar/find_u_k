install.packages("primes", repos = "http://cran.us.r-project.org")
install.packages("rlist", repos = "http://cran.us.r-project.org")
library(primes)
library(MASS)
library(rlist)

# Initial setup
primes <- primes::generate_primes(smallest_prime_to_consider, max_prime)
p_1 <- primes[1]

# Traverses the elements l in U_{k-1} to add l+1 to U_k, and 
# k-1+diff to U_k, where diff=k-l.
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
# cur_sum: current factorization sum.
# fact_length: current factorization length.
add_new_lengths_for_n_rec <- function(n, max_atom_index_to_consider, cur_sum, fact_length){
  if(n == cur_sum) return (c(fact_length))
  if(cur_sum > n) return (c(NA))
  result <- list()
  for(i in 1:max_atom_index_to_consider){
    temp <- add_new_lengths_for_n_rec(n, max_atom_index_to_consider, cur_sum + primes[i] - 1, fact_length + primes[i])
    result[[length(result) + 1]] <- temp
  }
  return (result)
}

# Adds the lengths not found with method add_known_lengths_from_u_k_minus_1.
# If $x$ has a factorization of orden $k$ not obtained from $U_{k-1}$,
# then $x$ must be an integer in certain interval. So it sufficies to iterate through such
# interval finding all the possible lenghts of elements in that interval.
# This method iterates through all the possible values of $x$.
add_new_lengths_to_u_k <- function(k){
  #k <- 5
  #n <- 4
  u_k <- list()
  max_integer <- k - 1
  min_integer <- ceiling(k*(p_1 - 1)/p_1)
  # Each atom must be multiplied by a multiple of its denominator.
  for(n in min_integer:max_integer){
    max_atom_index_to_consider <- max_atom_index_to_consider(n)
    if(max_atom_index_to_consider == 0 || (n < (primes[1] - 1))){# Do some basic checks
      print("Either max_atom_index_to_consider is 0 or n < (primes[1] - 1)")
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

# This method iterates through all possible lengths of u_k
# that are not in 'known_lengths_from_u_k_minus_1' and checks them
# against the 'new_lengths_in_u_k' vector.
# It returns TRUE if 'new_lengths_in_u_k' vector' is correct for the given k.
test_new_lengths_in_u_k <- function(k, known_lengths_from_u_k_minus_1, new_lengths_in_u_k){
  
  # Testing
  #k <- 3
  #known_lengths_from_u_k_minus_1 <- c(3)
  #new_lengths_in_u_k <- c(3, 4)
  
  lower_bound_u_k <- ceiling(k*(p_1 - 1)/p_1)
  upper_bound_u_k <- floor(k*p_1/(p_1 - 1))
  tentative_u_k <- lower_bound_u_k:upper_bound_u_k
  
  for(l in tentative_u_k){
    
    # Testing
    #l <- 2
    
    if(!l %in% known_lengths_from_u_k_minus_1 && l != k){
      indexes <- is_l_in_u_k(k, l)
      if(length(indexes) == 0 && l %in% new_lengths_in_u_k)
      {
        print(paste0("The following length should NOT be in u_", k, ": ", l))
        return(FALSE)
      }
      if(length(indexes) > 0 && !(l %in% new_lengths_in_u_k)){
        print(paste0("The following length is MISSING from u_", k, ": ", l))
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}

# Returns an empty vector if n doesn't have a factorization of length l.
# Otherwise, it returns the indexes of the atoms corresponding 
# to a factorization of n of lenght l. 
is_l_in_set_of_lengths_of_n_rec <- function(n, l, cur_sum, indexes){
  
  # Testing
  #k <- 3
  #l <- 4
  #n <- 2
  #cur_sum <- 0
  #indexes <- c()
  
  if(cur_sum == n){
    temp_length <- sum(primes[indexes])
    if(temp_length == l)
      return(indexes)
    else
      return(c())
  }
  else {
    if(cur_sum > n)  return(c())
    else {
      max_atom <- max_atom_index_to_consider(n)
      for(i in 1:max_atom){
        
        # Testing
        #i <- 1
        
        temp <- is_l_in_set_of_lengths_of_n_rec(n, l, cur_sum + primes[i] - 1, c(indexes, i))
        if(length(temp) > 0) return(temp)
      }
    }
  }
  return(c())
}

# Returns a vector of indexes if there is an integer
# having a factorization of length l or 
# an empty vector otherwise.
is_l_in_u_k <- function(k, l){
  
  # Testing
  #k <- 3
  #l <- 4
  
  m <- min(k, l)
  for(n in ceiling(m*(p_1 - 1)/p_1):floor(m)){
    
    # Testing 
    #n <- 2
    
    indexes_k <- is_l_in_set_of_lengths_of_n_rec(n, k, 0, c())
    indexes_l <- is_l_in_set_of_lengths_of_n_rec(n, l, 0, c())
    if(length(indexes_k) > 0 && length(indexes_l) > 0) return(indexes_l)
  }
  
  return(c())
}

# Finds u_k given u_{k-1}
find_u_k <- function(k, u_k_minus_1){
  
  # Testing
  #k <- 5
  #u_k_minus_1 <- c(3, 4, 5)
  
  
  known_lengths_from_u_k_minus_1 <- add_known_lengths_from_u_k_minus_1(k, u_k_minus_1)
  new_lengths_in_u_k <- add_new_lengths_to_u_k(k)
  
  #test <- test_new_lengths_in_u_k(k, known_lengths_from_u_k_minus_1, new_lengths_in_u_k)
  #if(!test)
  #{
  #  print("The values of u_k are not correct.")
  #  return(NA)
  #}
  
  u_k <- unique(unlist(c(known_lengths_from_u_k_minus_1, new_lengths_in_u_k)))
  
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
    print(paste("Done with k = ", sep = "", i))
  }
  names(set_of_unions) <- paste0("U_", 1:max_k)
  return(set_of_unions)
}

find_set_of_union_of_sets_of_lengths_from_min_k_to_max_k <- function(min_k, max_k, set_of_union_of_sets_of_lengths_for_k_less_than_min_k){
  #max_k = 3
  for(i in min_k:max_k){
    set_of_union_of_sets_of_lengths_for_k_less_than_min_k[[i]] <- find_u_k(i, set_of_union_of_sets_of_lengths_for_k_less_than_min_k[[i-1]])
  }
  names(set_of_union_of_sets_of_lengths_for_k_less_than_min_k) <- paste0("U_", 1:max_k)
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
       pch=18, col=2,  xlab="k", ylab="U_k", axes = FALSE)
  #plot(xValues, yValues, cex=dotsize, xlim=c(0,width+1), ylim=c(0, height+1), 
  #     bty="n", pch=1,  col=cm.colors(2), xlab="k", ylab="Uk", axes = FALSE)
  axis(side = 1, at = 1:width)
  axis(side = 2, at = 1:height)
}