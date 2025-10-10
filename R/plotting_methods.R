# Functions to plot the copula data and the parameter functions

#' Marginally normalized contour plots for copulas with different
#' families and one-dimensional conditioning values
#' @param cond_vals The conditioning values
#' (plotted next to each other in columns)
#' @param family_vals The names of the copula families in the rows
#' @param param_cond_func_1d The parameter function to use
#' @param automatic_title Whether to automatically create a title (TRUE or FALSE)
#' @param manual_title Manually enter a title if desired. If this is not "",
#' then automatic_title is automatically set to FALSE
#' @param manual_subtitle Manually enter a subtitle if desired. If this is "",
#' then an automatic subtitle is set,in
#' if either automatic_title is TRUE or manual_title is not "".
plot_contours_1d <- function(
    cond_vals=c(0.2,0.4,0.6,0.8),
    family_vals=c("frank","gaussian","gumbel", "clayton"),
    param_cond_func_1d = u_to_param_linear(c(1)),
    automatic_title=TRUE,
    manual_title="",
    manual_subtitle=""){
  n <- length(cond_vals)
  m <- length(family_vals)
  # leave enough space for a title on top of the plots, if desired
  space_on_top <- 1.5
  if(automatic_title || manual_title != ""){
    space_on_top <- 4
  }
  par(mfrow=c(m, n), mar=c(1,1,1,1), oma=c(1.5,1.5,space_on_top,1.5))
  for(i in  1:m){
    for(j in 1:n){
      temp_param <- param_cond_func_1d(cond_vals[j], family=family_vals[i])
      bicop_cont_dist <- bicop_dist(
        family=family_vals[i],
        rotation=0,
        parameters = temp_param)
      contour(bicop_cont_dist,margins="norm",
              main=paste(family_vals[i], round(cond_vals[j],2)),
              drawlabels=FALSE,
              axes=FALSE)
      box()
    }
  }
  if(automatic_title && manual_title == ""){
    mtext("Marginally normalized contour plots", outer = TRUE, cex = 1, line = 2, font=2)
    if(manual_subtitle==""){
      mtext("Different families and conditioning values, see titles", outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle, outer = TRUE, cex = 0.8, line = 1)
    }
  }
  if(manual_title != ""){
    mtext(manual_title, outer=TRUE, cex=1, line=2, font=2)
    if(manual_subtitle==""){
      mtext("Different families and conditioning values, see titles", outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle, outer = TRUE, cex = 0.8, line = 1)
    }
  }
}

#' Marginally normalized contour plots for a copula with family family_name
#' with different two-dimensional conditioning values
#' @param cond_vals_1 The conditioning values of the first variable
#' (plotted next to each other in columns)
#' @param cond_vals_2 The conditioning values of the second variable
#' (plotted in the rows)
#' @param family_name The name of the copula family to use
#' @param param_cond_func_2d The parameter function to use
#' @param automatic_title Whether an automatically created title
#' should be plotted.Needs to be TRUE or FALSE.
#' @param manual_title Manual Title, if desired. If this is not "", then
#' automatic_title is automatically set to FALSE.
#' @param manual_subtitle Manual Subtitle, if desired. If this is "",
#' then an automatic subtitle is set,
#' if either automatic_title is TRUE or manual_title is not "".
plot_contours_2d <- function(
    cond_vals_1 = c(0.2,0.4,0.6,0.8),
    cond_vals_2 = c(0.2,0.4,0.6,0.8),
    family_name="frank",
    param_cond_func_2d = u_to_param_linear(c(0.7,0.3)),
    automatic_title = TRUE,
    manual_title="",
    manual_subtitle=""){
  m <- length(cond_vals_1)
  n <- length(cond_vals_2)
  # leave enough space for a title on top of the plots, if desired
  space_on_top <- 1.5
  if(automatic_title || manual_title != ""){
    space_on_top <- 4
  }
  par(mfrow=c(m, n), mar=c(1,1,1,1), oma=c(1.5,1.5,space_on_top,1.5))
  for(i in  1:m){
    for(j in 1:n){
      temp_param <- param_cond_func_2d(c(cond_vals_1[i], cond_vals_2[j]), family=family_name)
      bicop_cont_dist <- bicop_dist(
        family=family_name,
        rotation=0,
        parameters = temp_param)
      contour(bicop_cont_dist,margins="norm",
              main=paste("1:", cond_vals_1[i], "2:", cond_vals_2[j]),
              drawlabels=FALSE,
              axes=FALSE)
      box()
    }
  }
  if(automatic_title && manual_title == ""){
    mtext(paste0("Marginally normalized contour plots for the ", family_name, " copula"),
          outer = TRUE, cex = 1, line = 2, font=2)
    if(manual_subtitle==""){
      mtext("Two conditioning values, see plot titles",
            outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle,
            outer = TRUE, cex = 0.8, line = 1)
    }
  }
  if(manual_title != ""){
    mtext(manual_title, outer=TRUE, cex=1, line=2, font=2)
    if(manual_subtitle==""){
      mtext("Two conditioning values, see plot titles",
            outer = TRUE, cex = 0.8, line = 1)
    } else {
      mtext(manual_subtitle,
            outer = TRUE, cex = 0.8, line = 1)
    }
  }
}

#' Function with similar functionality as the pairs_copula_data function of
#' the rvinecopulib package. Plots histograms of data on the diagonal,
#' scatter plots on the upper triangle and
#' marginally normalized contours on the lower triangle.
#' @param data Data to plot
#' @param uscale True, if the data is already on the copula scale, False otherwise.
#' @param method Method to use for transforming the data to copula scale, if uscale is False.
#' @param title Title for the Plot.
#' @return An object, which can be plotted using ggplotify::as.ggplot(object)
#' @export
copula_pairs_ggplot <- function(data, uscale = TRUE, method = c("ecdf", "kde1d"), title = NULL) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(all(sapply(data, is.numeric)))
  stopifnot(!anyNA(data))
  data <- as.data.frame(data)
  method <- match.arg(method)
  lapply(c("ggplot2", "gridExtra", "grid", "dplyr"), requireNamespace, quietly = TRUE)

  var_names <- colnames(data)
  n <- ncol(data)
  text_size <- 8 / sqrt(n)

  if (uscale) {
    copula_data <- data
  } else {
    # --- Transform to copula scale ---
    copula_data <- switch(method,
                          ecdf = as.data.frame(apply(data, 2, function(x) rank(x) / (length(x) + 1))),
                          kde1d = {
                            if (!requireNamespace("kde1d", quietly = TRUE)) {
                              stop("Package 'kde1d' is required for method = 'kde1d'")
                            }
                            as.data.frame(data) |>
                              dplyr::mutate(dplyr::across(dplyr::everything(), ~ {
                                fit <- kde1d::kde1d(.x)
                                u <- kde1d::pkde1d(.x, fit)
                                pmin(pmax(u, 1e-6), 1 - 1e-6)
                              }))
                          }
    )
  }

  # Gaussian margins
  gaussian_data <- as.data.frame(lapply(copula_data, qnorm))
  colnames(gaussian_data) <- var_names

  # --- Diagonal histogram plot ---
  plot_hist_diag <- function(i) {
    var_data <- copula_data[[var_names[i]]]
    hist_data <- hist(var_data, breaks = seq(0, 1, length.out = 31), plot = FALSE)
    y_max <- max(hist_data$counts) * 1.5

    ggplot2::ggplot(copula_data, ggplot2::aes(x = !!rlang::sym(var_names[i]))) +
      ggplot2::geom_histogram(breaks = seq(0, 1, length.out = 31),
                              fill = "grey80", color = "black") +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, y_max), expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::annotate("text", x = 0.5, y = y_max, label = var_names[i],
                        hjust = 0.5, vjust = 1, size = text_size) +
      ggplot2::theme_void()
  }

  # --- Upper triangle: scatter/hex + tau ---
  plot_scatter_upper <- function(i, j) {
    tau <- cor(copula_data[[j]], copula_data[[i]], method = "kendall")
    ggplot2::ggplot(copula_data, ggplot2::aes(x = !!rlang::sym(var_names[j]), y = !!rlang::sym(var_names[i]))) +
      ggplot2::geom_hex(alpha = 0.4, bins = 30) +
      ggplot2::annotate("text", x = 0.05, y = 0.95,
                        # label = format(round(tau, 2), digits = 2),
                        label = sprintf("%.2f", tau),
                        size = text_size, hjust = 0, vjust = 1, color = "red") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")
  }

  # --- Lower triangle: contour ---
  plot_contour_lower <- function(i, j) {
    x <- gaussian_data[[j]]
    y <- gaussian_data[[i]]
    valid <- all(is.finite(x)) && all(is.finite(y)) &&
      length(unique(x)) > 5 && length(unique(y)) > 5 &&
      sd(x) > 1e-8 && sd(y) > 1e-8

    if (!valid) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No variation",
                          size = 4, hjust = 0.5) +
        ggplot2::theme_void()
    } else {
      ggplot2::ggplot(data.frame(x = x, y = y), ggplot2::aes(x = x, y = y)) +
        ggplot2::stat_density_2d(
          ggplot2::aes(fill = ggplot2::after_stat(level)),
          geom = "polygon", bins = 5, color = NA
        ) +
        ggplot2::scale_fill_viridis_c(option = "D") +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none")
    }
  }

  # --- Build matrix of plots ---
  plot_matrix <- vector("list", n * n)
  for (i in 1:n) {
    for (j in 1:n) {
      idx <- (i - 1) * n + j
      plot_matrix[[idx]] <- if (i == j) {
        plot_hist_diag(i)
      } else if (i < j) {
        plot_scatter_upper(i, j)
      } else {
        plot_contour_lower(i, j)
      }
    }
  }

  g_body <- gridExtra::arrangeGrob(grobs = plot_matrix, ncol = n)

  g_final <- if (!is.null(title)) {
    gridExtra::arrangeGrob(g_body, top = grid::textGrob(title, gp = grid::gpar(fontsize = 14, fontface = "bold")))
  } else {
    g_body
  }
  #grid::grid.draw(g_final)
  return(g_final)
}

#' Takes a vector a and returns a function of u (vector of elements between 0 and 1).
#' That function takes the dot product between a and u, scales that to
#' T_values and calculates a kendall's tau value using the tanh.
#' @param a a vector which needs to be a convex combination
#' (i.e. all entries >=0 and they have to sum to 1)
#' @param tau_lower Lowest tau value, defaults to -0.92
#' (should be between -1 and 1 and less than tau_upper)
#' @param tau_upper Highest tau value, defaults to 0.92
#' (should be between -1 and 1 and greater than tau_lower)
#' @return A function of u, which calculates the kendalls tau value of
#' a copula given the conditioned values.
u_to_ktau_linear <- function(a, tau_lower=-0.92, tau_upper=0.92){
  return(function(u){
    tryCatch({
      T_upper <- fisher_z_transform(tau_upper)
      T_lower <- fisher_z_transform(tau_lower)
      T_val <- T_lower
      if(length(u) == 1){
        arg <- u
        T_val <- (T_upper- T_lower) * arg + T_lower
      } else {
        arg <- a %*% u
        T_val <- (T_upper - T_lower) * arg + T_lower
      }
      tau <- inverse_fisher_transform(T_val)
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
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
u_to_ktau_quadratic <- function(a, tau_lower=-0.92, tau_upper=0.92){
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
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
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
u_to_ktau_cubic <- function(a, tau_lower=-0.92, tau_upper=0.92){
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
      return(tau)
    },
    error = function(e){
      stop(paste0("An error occurred:", e))
    })
  })
}
#' Plot of the cond to ktau function with 1 conditioning variable
#' @param func_list list of functions for which plots should be created
#' @param titles vector of titles for the plots corresponding to the different
#' functions in func_list
#' @param u_cond_vals Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values, for which the functions should be evaluated for the plot.
plot_cond_to_ktau_1d <- function(func_list, titles=-1, u_cond_vals=1:99/100){
  par(mfrow=c(length(func_list), 1),
      mar = c(3, 3, 1.5, 1.5),
      oma = c(1.5, 0, 1.5, 0),
      mgp=c(2,1,0))
  for(i in 1:length(func_list)){
    ktau_vals <- sapply(u_cond_vals, func_list[[i]])
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    plot(u_cond_vals,
         ktau_vals,
         type="l",
         main=plot_title,
         xlab="u",
         ylab="ktau",
         cex.lab=0.9)
  }
  par(mfrow=c(1,1))
}

#' 3d plot of the cond to ktau function, for 2 conditioning variables
#' @param func_list list of functions for which plots should be created
#' @param titles vector of titles for the plots corresponding to the different
#' functions in func_list
#' @param u_cond_vals_1 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the first conditioning variable, for which the
#' functions should be evaluated for the plot.
#' @param u_cond_vals_2 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the second conditioning variable, for which the
#' functions should be evaluated for the plot.
plot_cond_to_ktau_2d <- function(func_list,
                                 titles=-1,
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100){
  par(mfrow=c(1, length(func_list)), mar=c(1,1,1,1), oma=c(1.5,1.5,4,1.5))
  for(i in 1:length(func_list)){
    ktau_vals <- matrix(
      rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
      ncol=length(u_cond_vals_2))
    for(j in 1:length(u_cond_vals_1)){
      for(k in 1:length(u_cond_vals_2)){
        ktau_vals[j,k] <- func_list[[i]](c(u_cond_vals_1[j], u_cond_vals_2[k]))
      }
    }
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    # Create a 3D surface plot
    blue_gradient <- colorRampPalette(c("lightblue", "darkblue"))
    plot3D::persp3D(u_cond_vals_1, u_cond_vals_2, ktau_vals,
                    theta = 30, phi = 20, axes=TRUE, ticktype="detailed",
                    xlab="u1", ylab="u2", zlab="ktau", colvar=ktau_vals,
                    col = blue_gradient(100), border = "grey", lwd=0.1, main=plot_title)
  }
}

#' Contours of the cond to ktau function with 2 conditional variables
#' @param func_list list of functions for which plots should be created
#' @param titles vector of titles for the plots corresponding to the different
#' functions in func_list
#' @param u_cond_vals_1 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the first conditioning variable, for which the
#' functions should be evaluated for the plot.
#' @param u_cond_vals_2 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the second conditioning variable, for which the
#' functions should be evaluated for the plot.
plot_ktau_contours_2d <- function(func_list,
                                  titles=-1,
                                  u_cond_vals_1=1:99/100,
                                  u_cond_vals_2=1:99/100){
  for(i in 1:length(func_list)){
    ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                        ncol=length(u_cond_vals_2))
    for(j in 1:length(u_cond_vals_1)){
      for(k in 1:length(u_cond_vals_2)){
        ktau_vals[j,k] <- func_list[[i]](c(u_cond_vals_1[j], u_cond_vals_2[k]))
      }
    }
    plot_title <- "Unknown"
    if(all(titles != -1) & i <= length(titles)){
      plot_title <- titles[i]
    }
    # Create a contour plot
    contour(u_cond_vals_1, u_cond_vals_2, ktau_vals, main=plot_title,
            xlab="u1", ylab="u2")
  }
}


#' Several 3d Plots for in total 3 conditioning variables,
#' always one is fixed for the plots.
#' @param func_list list of functions for which plots should be created.
#' @param titles vector of titles for the plots corresponding to the different
#' functions in func_list.
#' @param u_cond_vals_1 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the first conditioning variable that is not fixed,
#' for which the functions should be evaluated for the plot.
#' @param u_cond_vals_2 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the second conditioning variable that is not fixed,
#' for which the functions should be evaluated for the plot.
#' @param u_cond_vals_fixed Vector, defaults to c(0.3,0.7). Contains the
#' conditioning values of the fixed variable.
plot_cond_to_ktau_3d <- function(cond_to_ktau_func,
                                 title="",
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100,
                                 u_cond_vals_fixed=c(0.25,0.5,0.75)){
  # order of mar and oma: bottom, left, top, right. oma=outer margin, mar = inner margin
  par(mfrow = c(3, length(u_cond_vals_fixed)),
      mar = c(1, 1, 1, 1),
      oma = c(1, 1, 4, 1))
  for(fixed_dim in 1:3){
    for(i in 1:length(u_cond_vals_fixed)){
      ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                          ncol=length(u_cond_vals_2))
      for(j in 1:length(u_cond_vals_1)){
        for(k in 1:length(u_cond_vals_2)){
          if(fixed_dim==1){
            # dimension 1 is fixed
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_fixed[i], u_cond_vals_1[j], u_cond_vals_2[k]))
          } else if(fixed_dim==2){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_fixed[i], u_cond_vals_2[k]))
          } else if(fixed_dim==3){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_2[k], u_cond_vals_fixed[i]))
          }
        }
      }
      # Determine title and axis labels
      plot_title <- paste0(
        "u", fixed_dim,"=", u_cond_vals_fixed[i])
      xlab_text <- "u1"
      ylab_text <- "u2"
      if(fixed_dim == 1){
        xlab_text <- "u2"
        ylab_text <- "u3"
      } else if(fixed_dim==2){
        xlab_text <- "u1"
        ylab_text <- "u3"
      }
      blue_gradient <- colorRampPalette(c("lightblue", "darkblue"))
      # Create a 3d plot
      plot3D::persp3D(u_cond_vals_1, u_cond_vals_2, ktau_vals,
                      theta = 30, phi = 20, axes=TRUE, ticktype="detailed",
                      cex.axis=0.5, cex.lab =0.7,
                      xlab=xlab_text, ylab=ylab_text,zlab="ktau",colvar=ktau_vals,
                      col = blue_gradient(100), border = "grey",lwd=0.1, main=plot_title)
    }
  }
  mtext(title, outer=TRUE, cex=1, line=2, font=2)
}


#' Contour plot of ktau for 3 conditioning variables, always one fixed.
#' @param func_list list of functions for which plots should be created.
#' @param titles vector of titles for the plots corresponding to the different
#' functions in func_list.
#' @param u_cond_vals_1 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the first conditioning variable that is not fixed,
#' for which the functions should be evaluated for the plot.
#' @param u_cond_vals_2 Vector, defaults to c(0.01,...,0.99). Contains the
#' conditioning values of the second conditioning variable that is not fixed,
#' for which the functions should be evaluated for the plot.
#' @param u_cond_vals_fixed Vector, defaults to c(0.25,0.5,0.75). Contains the
#' conditioning values of the fixed variable.
plot_ktau_contour_3d <- function(cond_to_ktau_func,
                                 title="",
                                 u_cond_vals_1=1:99/100,
                                 u_cond_vals_2=1:99/100,
                                 u_cond_vals_fixed=c(0.25,0.5,0.7)){
  par(mfrow = c(3, length(u_cond_vals_fixed)),
      mar = c(1.5, 1.5, 1.5, 1.5),
      oma = c(1.5, 4, 4, 4))
  for(fixed_dim in 1:3){
    for(i in 1:length(u_cond_vals_fixed)){
      ktau_vals <- matrix(rep(0,length(u_cond_vals_1)*length(u_cond_vals_2)),
                          ncol=length(u_cond_vals_2))
      for(j in 1:length(u_cond_vals_1)){
        for(k in 1:length(u_cond_vals_2)){
          if(fixed_dim==1){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_fixed[i], u_cond_vals_1[j], u_cond_vals_2[k]))
          } else if(fixed_dim==2){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_fixed[i], u_cond_vals_2[k]))
          } else if(fixed_dim==3){
            ktau_vals[j,k] <- cond_to_ktau_func(
              c(u_cond_vals_1[j], u_cond_vals_2[k], u_cond_vals_fixed[i]))
          }
        }
      }
      # Adjust title and axis labels
      plot_title <- paste0(
        "u", fixed_dim,"=", u_cond_vals_fixed[i])
      xlab_text <- "u1"
      ylab_text <- "u2"
      if(fixed_dim == 1){
        xlab_text <- "u2"
        ylab_text <- "u3"
      } else if(fixed_dim==2){
        xlab_text <- "u1"
        ylab_text <- "u3"
      }
      # Draw a contour plot
      contour(u_cond_vals_1, u_cond_vals_2, ktau_vals,
              xlab=xlab_text, ylab=ylab_text, main=plot_title)
    }
  }
  mtext(title, outer=TRUE, cex=1, line=2, font=2)
}
