# This script contains functions to simulate from a Non-Simplified Vine Copula

#' Computes the fisher Z transform (artanh)
#' @param tau the value of Kendall's tau to transform to a real number.
#' Can also be a vector of taus.
#' @return The fisher transformation(s) of tau. If tau is a number, this returns a number
#' if tau was a vector this returns a vector
fisher_z_transform = function(tau){
  return(atanh(tau))
}
#' Computes the inverse of the fisher Z transform, which is the tanh
#' @param T_val The T-value that is supposed to be turned into a Kendall's tau again
#' @return The inverse fisher transformation(s) of T_val. If T is a number, this returns a number
#' if T is a vector this returns a vector
inverse_fisher_transform <- function(T_val){
  # computes the inverse of the fisher Z transform (tanh)
  return(tanh(T_val))
}

#' Returns scaling_factor * tanh(T_val) + shift
#' @param T_val The T-value that is supposed to be turned into a Kendall's tau again
#' @param scaling_factor scale the result by this factor if desired. By default,
#' tanh returns values between -1 and 1. if a scaling_factor is specified, the result is
#' between - scaling factor and scaling factor
#' @param shift shift the result if desired
#' @return The inverse fisher transformation(s) of T_val. If T is a number, this returns a number
#' if T is a vector this returns a vector
scaled_tanh <- function(T_val, scaling_factor=1.0, shift=0.0){
  # computes the inverse of the fisher Z transform (tanh)
  return(scaling_factor * tanh(T_val) + shift)
}

#' Maps a conditioning value to a copula parameter. Can only map to two different values (
#' step function)
#' @param u the samples u value (between 0 and 1)
#' @param family the family of the third copula in the sample
#' @return A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
#' @export
u_to_param <- function(u, family="gaussian"){
  tryCatch({
    # tau should be in the range -0.92, +0.92 to avoid numerical problems later on
    tau_min <- -0.92
    tau_max <- 0.92
    T_min <- fisher_z_transform(tau_min)
    T_max <- fisher_z_transform(tau_max)
    #arg = ifelse(length(u)==1,u,max(u))
    arg=ifelse(length(u)==1,u,mean(u))
    #T_val <- (T_max - T_min)*arg + T_min
    T_val <- ifelse(arg > 0.5,T_max,T_min)
    # T_val <- 4*(T_max - T_min) * ((arg-0.5)^2) + T_min
    tau <- inverse_fisher_transform(T_val)
    param <- ktau_to_par(family=family, tau=tau)
    return(param)
  },
  error = function(e){
    stop(pasteo0("An error occurred: ",
    e,
    ". A typical cause of errors:
    The function u_to_param only supports families for which ktau_to_par is defined.
    This error likely is thrown, because simulate_non_simplified was called with
    the default argument for param_cond_funcs with a family,
    for which ktau_to_par is not implemented."))
  })
}

#' takes two floats tau_lower and tau_upper and returns a function of u
#' (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function is constant and always returns 0.5 * (tau_lower + tau_upper)
#' @param tau_lower Lowest tau value, defaults to -0.92
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper Highest tau value, defaults to 0.92
#' (should be between -1 and 1 and greater than tau_lower)
#' @return A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
#' @export
u_to_param_constant <- function(tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      tau <- 0.5 * (tau_lower + tau_upper)
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        tau <- 0.5 * (tau_upper + 0.001)
      }
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          This function only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function takes the dot product between a and u
#' and returns a linear function of that dot product.
#' @param a a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
#' @param tau_lower Lowest tau value, defaults to -0.92
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper Highest tau value, defaults to 0.92
#' (should be between -1 and 1 and greater than tau_lower)
#' @return A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
#' @export
u_to_param_linear <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      T_upper <- fisher_z_transform(tau_upper)
      T_lower <- fisher_z_transform(tau_lower)
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        T_lower <- fisher_z_transform(0.001)
      }
      T_val <- T_lower
      if(length(u) == 1){
        arg <- u
        T_val <- (T_upper- T_lower) * arg + T_lower
      } else {
        arg <- a %*% u
        T_val <- (T_upper - T_lower) * arg + T_lower
      }
      tau <- inverse_fisher_transform(T_val)
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          This function only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function),
#' builds polynomial terms of degree 2 of these terms and then transforms them
#' to the interval (tau_lower, tau_upper) using a scaled tanh.
#' If dim(u) <=3 includes all terms u[i]*u[j], (powers and mixed terms).
#' If dim(u) > 3 only includes the individual powers, u[i]^k, k=1,2.
#' @param a a vector of weights
#' @param tau_lower Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @return A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
#' @export
u_to_param_quadratic <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      # define parameters for the scaled tanh, so that the result will be between
      # tau_lower and tau_upper
      scaling_factor <- (tau_upper - tau_lower) / 2
      shift <- (tau_upper + tau_lower) / 2
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        scaling_factor <- (tau_upper - 0.001) / 2
        shift <- (tau_upper + 0.001) / 2
      }
      arg <- qnorm(u)
      tau <- 0
      if(length(u) == 1){
        temp <- c(arg, arg^2) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==2) {
        temp <- c(arg[1], arg[2], arg[1]^2, arg[2]^2, arg[1]*arg[2]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==3){
        temp <- c(arg[1], arg[2], arg[3],
                  arg[1]^2, arg[2]^2, arg[3]^2,
                  arg[1]*arg[2], arg[1]*arg[3], arg[2]*arg[3]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      }else{
        # do not include mixed terms here.
        arg <- poly(arg, 2, raw=TRUE) # raw to just evaluate the polynomial terms
        # only keep the entries and list them in a vector with
        # c(arg[1],arg[2],..., arg[1]^2, arg[2]^2,...)
        arg <- c(unname(arg[1:nrow(arg),]))
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          This function only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1)
#' and of a string "family", which corresponds to a copula family.
#' that function transforms the values u using qnorm (inverse normal dist.function),
#' builds polynomial terms of degree 2 of these terms and then transforms them
#' to the interval (tau_lower, tau_upper) using a scaled tanh.
#' If dim(u) <=3 includes powers and mixed terms,
#' If dim(u) > 3 only includes the individual powers, u[i]^k, k=1,2,3.
#' @param a a vector of weights
#' @param tau_lower Lowest tau value
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper Highest tau value
#' (should be between -1 and 1 and greater than tau_lower)
#' @return A function of u and family, which calculates the parameter of
#' a copula given the conditioned values.
#' @export
u_to_param_cubic <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u, family="gaussian"){
    tryCatch({
      # define parameters for the scaled tanh, so that the result will be between
      # tau_lower and tau_upper
      scaling_factor <- (tau_upper - tau_lower) / 2
      shift <- (tau_upper + tau_lower) / 2
      if(family %in% c("clayton", "gumbel", "joe") && tau_lower <= 0.001){
        scaling_factor <- (tau_upper - 0.001) / 2
        shift <- (tau_upper + 0.001) / 2
      }
      arg <- qnorm(u)
      tau <- 0
      if(length(u) == 1){
        temp <- c(arg, arg^2, arg^3) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==2) {
        temp <- c(arg[1], arg[2],
                  arg[1]^2, arg[2]^2, arg[1]*arg[2],
                  arg[1]^3, arg[2]^3,arg[1]^2 * arg[2], arg[1]*arg[2]^2) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      } else if(length(u)==3){
        temp <- c(arg[1], arg[2], arg[3],
                  arg[1]^2, arg[2]^2, arg[3]^2,
                  arg[1]*arg[2], arg[1]*arg[3], arg[2]*arg[3],
                  arg[1]^3, arg[2]^3, arg[3]^3,
                  arg[1]^2*arg[2], arg[1]^2*arg[3],
                  arg[2]^2*arg[1], arg[2]^2*arg[3],
                  arg[3]^2*arg[1], arg[3]^2*arg[2], arg[1]*arg[2]*arg[3]) %*% a
        tau <- scaled_tanh(temp, scaling_factor=scaling_factor, shift=shift)
      }else{
        # do not include mixed terms here.
        arg <- poly(arg, 3, raw=TRUE) # raw to just evaluate the polynomial terms
        # only keep the entries and list them in a vector with
        # c(arg[1],arg[2],..., arg[1]^2, arg[2]^2,..., arg[1]^3,arg[2]^3,...)
        arg <- c(unname(arg[1:nrow(arg),]))
        tau <- scaled_tanh(a%*%arg, scaling_factor=scaling_Factor, shift=shift)
      }
      param <- ktau_to_par(family=family, tau=tau)
      return(param)
    },
    error = function(e){
      stop(paste0("An error occurred:", e, ". Common causes of an error:
          This function only supports families for which ktau_to_par is defined.
          The vector a does not have sufficiently many entries."))
    })
  })
}

#' permutes indices of a rvine structure matrix so that on the antidiagonal
#' the values 1 to d (dimension of the matrix) are ordered from top right to bottom left
#' @param mat The rvine structure matrix
#' @return list with matrix as first entry, another list as second entry and a vector as third entry
#' First entry of result is the new structure matrix, second works as a dictionary, mapping the values
#' the last contains the values on the antidiagonal of the original matrix.
permute_indices <- function(mat){
  antidiag_values <- rev(as_rvine_structure(mat)$order)
  n <- nrow(mat)
  permutation_dict <- list()
  for (i in 1:n){
    permutation_dict[[antidiag_values[i]]] <- i
  }
  struct_mat_new <- matrix(rep(0, length(mat)), nrow=nrow(mat))
  for (i in 1:n){
    for (j in 1:(n-i+1)){
      struct_mat_new[i,j] <- permutation_dict[[mat[i,j]]]
    }
  }
  return(list(struct_mat_new, permutation_dict, antidiag_values))
}

#' From the samples of the permuted structure, this function returns the samples from the original structure
#' @param simulations Matrix with ncol(simulations) = length(original_order)
#' The samples from the reordered rvine structure (with antidiagonal values ordered
#' from bottom left to top right)
#' @param original_order The elements on the antidiagonal from bottom left to top right in the original structure
#' @return result: the simulations, with the correct order.
permute_back <- function(simulations, original_order){
  result = simulations[,1:ncol(simulations)]
  for (i in 1:ncol(simulations)){
    result[,original_order[i]] <- simulations[,i]
  }
  return(result)
}

#' function to compute the matrix tilde(M)
#' @param mat The input matrix
#' @return the matrix, with entries tilde(m)_k,i = max_j=1,...,k (m_j,i)
#' @export
get_max_matrix <- function(mat){
  num_cols <- ncol(mat)
  num_rows <- nrow(mat)
  max_mat <- matrix(rep(0, num_cols * num_rows), ncol=num_cols)
  for (j in 1:num_cols){
    max_in_col <- -1
    for (i in 1:(num_rows -j+1)){
      if(mat[i,j] > max_in_col){
        max_in_col <- mat[i,j]
      }
      max_mat[i,j] <- max_in_col
    }
  }
  return(max_mat)
}

#' Simulates a single sample,
#' This function is called by simulate_non_simplified, after modifying the input structure,
#' to be able to parallelize the procedure more easily for better performance
#' @param new_struct rvine structure, where the antidiagonal is ordered
#' increasingly from top right to bottom left.
#' @param M_mat the matrix \eqn{\tilde{M}}, with \eqn{\tilde{m}_{k,i} = max_{j=1}^{k}{new_struct_{j,i}}}
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @return u_row: a d-dimensional vector sampled from the non-simplified vine copula
simulate_single_sample <- function(new_struct,
                                   M_mat,
                                   families,
                                   params,
                                   param_cond_funcs,
                                   rotations){
  d <- ncol(new_struct)
  w_row <- runif(d, min=0, max=1)
  u_row <- rep(0,d)
  # storage for the matrices V and V^2
  V <- matrix(rep(0,d^2), ncol=d)
  V_2 <- matrix(rep(0,d^2), ncol=d)
  V[1,d] <-w_row[d]
  # loop through the columns from right to left
  for (i in (d-1):1){
    V[d-i+1,i] <- w_row[i]
    # loop through the rows from the antidiagonal upwards
    for (k in (d-i):1){
      # In this iteration, we condition on m_{1,i},...,m_{k-1,i}
      # So for a non-simplified vine, the parameters can depend on
      # u_{m_{1,i}},..., u_{m_{k-1,i}}. These are already known at this point,
      # and they are stored in V_{1,{d-m_{1,i}+1}},...,V_{1,{d-m_{k-1,i}+1}}
      param <- params[[i]]
      if (k>1){
        cond_u_vals <- rep(0,k-1)
        for (temp in 1:(k-1)){
          cond_u_vals[temp] <- V[1,d-new_struct[temp,i]+1]
        }
        param <- param_cond_funcs[[k-1]][[i]](cond_u_vals, family=families[[k]][[i]])
      }
      cop <- bicop_dist(family=families[[k]][[i]],
                        rotation=rotations[[k]][[i]],
                        parameters = param)
      if(new_struct[k,i] == M_mat[k,i]){
        V[k,i] <- hbicop(c(V[k,d-M_mat[k,i]+1],
                           V[k+1,i]),
                         cop,
                         cond_var = 1, #i > m_{k,i} always holds for k <= d-i
                         inverse=TRUE)
      } else {
        V[k,i] <- hbicop(c(V_2[k, d-M_mat[k,i]+1],
                           V[k+1,i]),
                         cop,
                         cond_var = 1, #i > m_{k,i} always holds for k <= d-i
                         inverse=TRUE)
      }
      if(i>1){
        if(new_struct[k,i] == M_mat[k,i]){
          V_2[k+1,i] = hbicop(c(V[k,d-M_mat[k,i] + 1],
                                V[k,i]),
                              cop,
                              cond_var = 2, #i > m_{k,i} always holds for k <= d-i
                              inverse=FALSE)
        } else {
          V_2[k+1,i] = hbicop(c(V_2[k,d-M_mat[k,i] + 1],
                                V[k,i]),
                              cop,
                              cond_var=2, #i > m_{k,i} always holds for k <= d-i
                              inverse=FALSE)
        }
      }
    }
  }
  for(idx in 1:d){
    u_row[idx] <- V[1,d-idx+1]
  }
  return(u_row)
}

#' Simulation of a non_simplified vine copula.
#' For faster sampling use simulate_non_simp_parallel instead.
#' @param n_samples number of samples to simulate
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @return result: A nxd matrix, where each row contains one sample.
#' @export
simulate_non_simplified <- function(n_samples = 500,
                                    struct = matrix(c(1,1,1,
                                                      2,2,0,
                                                      3,0,0)
                                                    ,byrow=TRUE, ncol=3),
                                    families=list(list("gumbel", "gumbel"), list("gumbel")),
                                    params = list(c(2), c(1.3)),
                                    param_cond_funcs = list(list(u_to_param)),
                                    rotations = list(list(0,0),list(0))){
  # Reorder the indices, so they have the numbers 1 to d from top right to bottom left on the antidiagonal
  temp <- permute_indices(struct)
  new_struct <- temp[[1]]
  permutation <- temp[[2]]
  inverse_permutation <- temp[[3]]
  # Compute the matrix \tilde(M), with \tilde(M)_{k,i} = max(new_struct_{1,i},...,new_struct_{k,i})
  M_mat <- get_max_matrix(new_struct)
  # number of dimensions
  d <- ncol(new_struct)
  # matrix to save the simulations in
  u <- matrix(rep(0,d*n_samples),ncol=d)
  # loop once for every sample
  prog_bar <- txtProgressBar(min=0,max=n_samples,title="Simulating...", style=3)
  for (sample in 1:n_samples){
    u_sample <- simulate_single_sample(new_struct,
                           M_mat,
                           families,
                           params,
                           param_cond_funcs,
                           rotations)
    u[sample,] <- u_sample
    setTxtProgressBar(prog_bar, value=sample, title="Simulating...")
  }
  # permute the indices back to get a sample from the original rvine-structure
  result <- permute_back(u, inverse_permutation)
  return(result)
}



#' parallel simulation of a non_simplified vine copula.
#' @param n_samples number of samples to simulate
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @return result: A nxd matrix, where each row contains one sample.
#' @export
simulate_non_simp_parallel <- function(n_samples = 500,
                                       struct = matrix(c(1,1,1,
                                                         2,2,0,
                                                         3,0,0)
                                                       ,byrow=TRUE, ncol=3),
                                       families=list(list("gumbel", "gumbel"), list("gumbel")),
                                       params = list(c(2), c(1.3)),
                                       param_cond_funcs = list(list(u_to_param)),
                                       rotations = list(list(0,0),list(0))){
  # Reorder the indices, so they have the numbers 1 to d from top right to bottom left on the antidiagonal
  temp <- permute_indices(struct)
  new_struct <- temp[[1]]
  permutation <- temp[[2]]
  inverse_permutation <- temp[[3]]
  # Compute the matrix \tilde(M), with \tilde(M)_{k,i} = max(new_struct_{1,i},...,new_struct_{k,i})
  M_mat <- get_max_matrix(new_struct)
  # number of dimensions
  d <- ncol(new_struct)
  # matrix to save the simulations in
  u <- matrix(rep(0,d*n_samples),ncol=d)

  # parallel computing
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # export the required functions for parallel processing
  parallel::clusterExport(cl,
                varlist = c("simulate_single_sample",
                            "u_to_param",
                            "bicop_dist",
                            "hbicop",
                            "ktau_to_par",
                            "fisher_z_transform",
                            "inverse_fisher_transform",
                            "scaled_tanh",
                            "struct",
                            "families",
                            "params",
                            "param_cond_funcs",
                            "rotations"),
                envir=environment())
  # pblapply to display a progress bar
  output_list <- pbapply::pblapply(1:n_samples, function(i) {
    simulate_single_sample(new_struct,
                           M_mat,
                           families,
                           params,
                           param_cond_funcs,
                           rotations)
  }, cl=cl)
  # stop the cluster, to free resources on computer
  parallel::stopCluster(cl)

  u <- do.call(rbind, output_list)
  # permute the indices back to get a sample from the original rvine-structure
  result <- permute_back(u, inverse_permutation)
  return(result)
}

#' From data coming from an rvine structure with elements current_order on the top right to
#' bottom left of the antidiagonal, this function returns the samples of the permuted structure,
#' where the elements on the antidiagonal are ordered increasingly from top right to bottom left.
#' @param samples Matrix with ncol(simulations) = length(current_order)
#' The given data.
#' @param current_order The elements on the antidiagonal from bottom left to top right
#' in the original structure
#' @return result: the data permuted to fit an order, that is increasing from the top
#' right to the bottom left of the r vine structure matrix.
permute_forward <- function(samples, current_order){
  result = samples[,1:ncol(samples)]
  for (i in 1:ncol(samples)){
    result[,i] <- samples[,current_order[i]]
  }
  return(result)
}

#' Evaluates the likelihood of a single sample,
#' @param u_row Single sample (on copula scale)
#' @param new_struct rvine structure, where the antidiagonal is ordered
#' increasingly from top right to bottom left.
#' @param M_mat the matrix \eqn{\tilde{M}}, with \eqn{\tilde{m}_{k,i} = max_{j=1}^{k}{new_struct_{j,i}}}
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @return L: Log-Likelihood of the datapoint
log_likelihood_single_sample <- function(
    u_row,
    new_struct,
    M_mat,
    families,
    params,
    param_cond_funcs,
    rotations
    ){
  d <- ncol(new_struct)
  L <- 0
  # storage for the matrices V and V^2
  V <- matrix(rep(0,d^2), ncol=d)
  V_2 <- matrix(rep(0,d^2), ncol=d)
  for(i in 1:d){
    V[1,i] <- u_row[d-i+1]
  }
  # loop through the columns from right to left
  for (i in (d-1):1){
    # loop through the rows from the antidiagonal upwards
    for (k in 1:(d-i)){
      # In this iteration, we condition on m_{1,i},...,m_{k-1,i}
      # So for a non-simplified vine, the parameters can depend on
      # u_{m_{1,i}},..., u_{m_{k-1,i}}. These are already known at this point,
      # and they are stored in V_{1,{d-m_{1,i}+1}},...,V_{1,{d-m_{k-1,i}+1}}
      param <- params[[i]]
      if (k>1){
        cond_u_vals <- rep(0,k-1)
        for (temp in 1:(k-1)){
          cond_u_vals[temp] <- V[1,d-new_struct[temp,i]+1]
        }
        param <- param_cond_funcs[[k-1]][[i]](cond_u_vals, family=families[[k]][[i]])
      }
      cop <- bicop_dist(family=families[[k]][[i]],
                        rotation=rotations[[k]][[i]],
                        parameters = param)
      z_1 <- V[k,i]
      z_2 <- 0
      if(new_struct[k,i] == M_mat[k,i]){
        z_2 <- V[k,d-M_mat[k,i]+1]
      } else {
        z_2 <- V_2[k,d-M_mat[k,i]+1]
      }
      # z_1 = C_{d-i+1|m_{1,i},...,m_{k-1,i}}(...),
      # z_2 = C_{m_{k,i}|m_{1,i},...,m_{k-1,i}}(...),
      # so first z_2 then z_1 is the right order in the following copula,
      # since k <= d-i < d-i+1
      L <- L + log(dbicop(c(z_2,z_1) ,cop))
      V[k+1,i] = hbicop(c(z_1,z_2),cop,cond_var = 2, inverse=FALSE)
      V_2[k+1,i] <- hbicop(c(z_1,z_2), cop, cond_var=1, inverse=FALSE)
    }
  }
  return(L)
}


#' Calculate likelihood
#' @param u_data matrix of datapoints. Each row should by one datapoint, for which
#' a likelihood should be evaluated
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @param return_vector Boolean, defaults to FALSE. If TRUE, for every row
#' in u_data the function returns the likelihood for that row, i.e. returns a vector. Otherwise
#' it returns a single number (the sum of the individual log likelihoods)
#' @return L: if return_vector is TRUE a vector, otherwise a number.
#' The log likelihoods of the individual samples, or of the whole dataset
#' @export
log_likelihood_non_simplified <- function(
    u_data =matrix(c(0.2,0.3,0.4,
                     0.5,0.5,0.5),ncol=3, byrow=TRUE),
    struct = matrix(c(1,1,1,
                      2,2,0,
                      3,0,0)
                    ,byrow=TRUE, ncol=3),
    families=list(list("gumbel", "gumbel"), list("gumbel")),
    params = list(c(2), c(1.3)),
    param_cond_funcs = list(list(u_to_param)),
    rotations = list(list(0,0),list(0)),
    return_vector=FALSE){
  # Reorder the indices, so they have the numbers 1 to d from top right to bottom left on the antidiagonal
  temp <- permute_indices(struct)
  new_struct <- temp[[1]]
  current_order <- temp[[3]]
  u_data_ordered <- permute_forward(samples=u_data, current_order=current_order)
  # Compute the matrix \tilde(M), with \tilde(M)_{k,i} = max(new_struct_{1,i},...,new_struct_{k,i})
  M_mat <- get_max_matrix(new_struct)
  # save log likelihood
  L <-0
  if(return_vector){
    L <- rep(0, nrow(u_data_ordered))
  }
  prog_bar <- txtProgressBar(min=0,max=nrow(u_data_ordered),title="Evaluating...", style=3)
  for (row_idx in 1:nrow(u_data_ordered)){
    log_lik_current_sample <- log_likelihood_single_sample(
      u_row= u_data_ordered[row_idx,],
      new_struct = new_struct,
      M_mat = M_mat,
      families = families,
      params = params,
      param_cond_funcs = param_cond_funcs,
      rotations = rotations
    )
    setTxtProgressBar(prog_bar, value=row_idx, title="Evaluating...")
    if(return_vector){
      L[row_idx] <- log_lik_current_sample
    } else {
      L <- L + log_lik_current_sample
    }
  }
  return(L)
}


#' parallel evaluation of log likelihood data
#' @param u_data matrix of datapoints. Each row should by one datapoint, for which
#' a likelihood should be evaluated
#' @param struct rvine-structure matrix
#' @param families copula families, formatted as a list of lists
#' @param params parameters for the first tree, formatted as a list (of vectors)
#' @param param_cond_funcs list of list of functions that maps the u values, which the respective
#' copula is conditioned on, and a copula family-name to a parameter for that copula
#' @param rotations rotations of the copulas, formatted as a list of lists.
#' @param return_vector Boolean, defaults to FALSE. If TRUE, for every row
#' in u_data the function returns the likelihood for that row, i.e. returns a vector. Otherwise
#' it returns a single number (the sum of the individual log likelihoods)
#' @return L: if return_vector is TRUE a vector, otherwise a number.
#' The log likelihoods of the individual samples, or of the whole dataset
#' @export
log_likelihood_non_simp_parallel <- function(
    u_data = matrix(c(0.2,0.3,0.4,
                      0.5,0.5,0.5),ncol=3, byrow=TRUE),
    struct = matrix(c(1,1,1,
                      2,2,0,
                      3,0,0)
                    ,byrow=TRUE, ncol=3),
    families=list(list("gumbel", "gumbel"), list("gumbel")),
    params = list(c(2), c(1.3)),
    param_cond_funcs = list(list(u_to_param)),
    rotations = list(list(0,0),list(0)),
    return_vector=FALSE){
  # Reorder the indices, so they have the numbers 1 to d from top right to bottom left on the antidiagonal
  temp <- permute_indices(struct)
  new_struct <- temp[[1]]
  current_order <- temp[[3]]
  u_data_ordered <- permute_forward(samples=u_data, current_order=current_order)
  # Compute the matrix \tilde(M), with \tilde(M)_{k,i} = max(new_struct_{1,i},...,new_struct_{k,i})
  M_mat <- get_max_matrix(new_struct)
  # save log likelihood
  L <-0
  if(return_vector){
    L <- rep(0, nrow(u_data_ordered))
  }
  # parallel computing
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # export the required functions for parallel processing
  parallel::clusterExport(cl,
                          varlist = c("log_likelihood_single_sample",
                                      "u_to_param",
                                      "bicop_dist",
                                      "hbicop",
                                      "ktau_to_par",
                                      "fisher_z_transform",
                                      "inverse_fisher_transform",
                                      "scaled_tanh",
                                      "struct",
                                      "families",
                                      "params",
                                      "param_cond_funcs",
                                      "rotations"),
                          envir=environment())
  # pblapply to display a progress bar
  output_list <- pbapply::pblapply(1:nrow(u_data_ordered), function(i) {
    log_likelihood_single_sample(
      u_row=u_data_ordered[i,],
      new_struct= new_struct,
      M_mat=M_mat,
      families = families,
      params = params,
      param_cond_funcs = param_cond_funcs,
      rotations = rotations
      )
  }, cl=cl)
  # stop the cluster, to free resources on computer
  parallel::stopCluster(cl)
  if(return_vector){
    return(unlist(output_list))
  }
  L <- do.call(sum, output_list)
  return(L)
}

