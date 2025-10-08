# This script contains methods for classification and quantile regression.

#' Performs a train test split for the original data
#' @param data The data set
#' @param train_perc Defaults to 0.8. Percentage to use for training.
#' @return A list with 2 entries: orig_train, orig_test
#' @export
train_test_split_orig <- function(data, train_perc){
  num_train <- floor(train_perc * nrow(data))
  train_idx <- sample(nrow(data),num_train)
  data_train <- data[train_idx,]
  data_test <- data[-train_idx,]
  return(list(data_train, data_test))
}

#' Performs a train test split with original data already split in train and test
#' Gives all datapoints in orig_data a label of 1, all datapoints in simplified data a label of 0
#' @param orig_data_train The training set of the real data, which gets labels of 1
#' @param orig_data_test The test set of the real data, which gets labels of 1
#' @param simplified_data Synthetic data set, simulated from a simplified
#' vine copula, gets labels of 0
#' @return A list with 4 entries, x_train, x_test, y_train, y_test in that order.
#' @export
train_test_split <- function(orig_data_train, orig_data_test, simplified_data){
  train_perc <- nrow(orig_data_train) / (nrow(orig_data_train) + nrow(orig_data_test))
  y_train_orig <- as.matrix(rep(1L, nrow(orig_data_train)))
  y_test_orig <- as.matrix(rep(1L, nrow(orig_data_test)))
  simplified_labels <- as.matrix(rep(0L, nrow(simplified_data)))
  num_train_simp <- floor(train_perc*nrow(simplified_data))
  train_simp_idx <- sample(nrow(simplified_data), num_train_simp)
  x_train_simp <- simplified_data[train_simp_idx,]
  y_train_simp <- as.matrix(simplified_labels[train_simp_idx,])
  x_test_simp <- simplified_data[-train_simp_idx,]
  y_test_simp <- as.matrix(simplified_labels[-train_simp_idx,])
  # combine the data into one dataset and shuffle
  x_train <- rbind(orig_data_train, x_train_simp)
  y_train <- rbind(y_train_orig, y_train_simp)
  shuffle_idx_train <- sample(nrow(x_train), nrow(x_train))
  x_train <- x_train[shuffle_idx_train,]
  y_train <- y_train[shuffle_idx_train,]
  x_test <- rbind(orig_data_test, x_test_simp)
  y_test <- rbind(y_test_orig, y_test_simp)
  shuffle_idx_test <- sample(nrow(x_test), nrow(x_test))
  x_test <- x_test[shuffle_idx_test,]
  y_test <- y_test[shuffle_idx_test,]
  return(list(x_train,
              x_test,
              y_train,
              y_test))
}

#' learning rate scheduler: halves the learning rate every 30 epochs
#' argument to a keras function
#' @param epoch current epoch
#' @param lr current learning rate
#' @return lr: the updated learning rate
#' @export
lr_schedule_fun <- function(epoch, lr) {
  if((epoch + 1) %% 30 == 0){
    return(lr/2)
  } else {
    return(lr)
  }
}

#' Defines and compiles a neural network for binary classification.
#' Activation function is leaky_relu, output is one layer with no activation
#' (to train a binary classifier with from_logits=TRUE option)
#' @param input_dim dimension of the input, defaults to 5
#' @param hidden_units vector containing the number of units per layer.
#' Defaults to c(20,10)
#' @param initial_lr initial learing rate to use, defaults to 0.01
#' (íf no lr_scheduler is used during training, this stays the same during
#' the whole training process)
#' @param use_tanh If TRUE uses tanh activation, otherwise uses leaky_relu.
#' @param leaky_relu_alpha slope in the negative part of the relu. Defaults to 0.01.
#' @return the compiled neural network model
#' @export
build_model <- function(
    input_dim=5,
    hidden_units=c(20,10),
    initial_lr = 0.01,
    use_tanh=FALSE,
    leaky_relu_alpha=0.1){
  model <- -1 # reset the model variable.
  # definition of the model
  if(use_tanh){
    model <- keras_model_sequential() %>%
      layer_dense(units = hidden_units[1], input_shape = input_dim, activation='tanh') %>%
      layer_dense(units = hidden_units[2], activation='tanh') %>%
      layer_dense(units = 1)
  } else {
    model <- keras_model_sequential() %>%
      layer_dense(units = hidden_units[1], input_shape = input_dim) %>%
      layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
      layer_dense(units = hidden_units[2]) %>%
      layer_activation_leaky_relu(alpha = leaky_relu_alpha) %>%
      layer_dense(units = 1)
  }
  model %>% compile(
    optimizer = keras$optimizers$Adam(learning_rate=initial_lr),
    loss = keras$losses$BinaryCrossentropy(from_logits=TRUE),
    metrics = list(binary_accuracy_from_logits)
  )
  return(model)
}

#' function to train the model
#' @param model The NN model to train
#' @param x_train training data
#' @param y_train labels for the training data
#' @param lr_schedule Function for learning rate scheduling. Defaults to lr_schedule_fun.
#' @param num_epochs Integer, defaults to 500.
#' Number of epochs for which the model should be trained.
#' @param batch_size Integer, defaults to 128.
#' Batch size to use to compute gradients during training
#' @param val_split Number between 0 and 1. Fraction of data used for validation.
#' @param verbose Whether output should be printed. Can be 0,1 or 2.
#' @return list(model,history): a list of the trained model in the first position
#' and the training history (loss, accuracy and learing_rate per epoch)
#' @export
train_model <- function(
    model,
    x_train,
    y_train,
    lr_schedule = lr_schedule_fun,
    num_epochs = 300,
    batch_size=128,
    val_split = 0.2,
    verbose=1
){
  # create the keras object required for learning rate scheduling
  lr_scheduler <- callback_learning_rate_scheduler(schedule = lr_schedule)

  # train the model
  history <- model %>% fit(
    x_train, y_train,
    epochs = num_epochs,
    batch_size = batch_size,
    validation_split=val_split, # use 20 percent of training data as validation data
    verbose = verbose, # 0 for slightly faster training (no output), 1 to observe progress while training
    callbacks=list(lr_scheduler)
  )
  return(list(model, history))
}

# METHODS FOR QUANTILE REGRESSION

#' Calculate the correction factors, i.e. the factors lambda, with which the
#' simplified copula needs to be multiplied to get the NCE copula.
#' @param model The classifier trained to distinguish non-simplified data
#' with labels 1 from simplified data with labels 0
#' @param obs the observations, for which the factors should be created.
#' @param nu The fraction of noise to real samples
#' @param verbose Whether output should be printed. Can be 0,1 or 2.
#' @return The correction factors
#' @export
correction_factors <- function(model, obs, nu=1, verbose=1){
  predictions <- model %>% predict(obs, verbose=verbose)
  # cast outputs of neural net to numeric, just to be sure
  cor_factors <- exp(as.numeric(predictions) + log(nu))
  return(cor_factors)
}

#' compute the value of the NCE copula obtained from the simplified fit,
#' together with the classifier
#' @param model The classifier trained to distinguish non-simplified data
#' with labels 1 from simplified data with labels 0
#' @param fitted_vine The vine copula fitted to the observed data.
#' @param obs The observations, for which the NCE copula should be created.
#' @param nu The fraction of noise to real samples
#' @param verbose Whether output should be printed. Can be 0,1 or 2.
#' @return The likelihoods of the copula derived via noise contrastive estimation
#' evaluated at the points given by obs.
#' @export
NCE_cop <- function(model, fitted_vine, obs, nu=1, verbose=1){
  cor_facs <- correction_factors(model=model, obs=obs, nu=nu, verbose=verbose)
  c_nce <- cor_facs * dvinecop(obs, fitted_vine)
  return(c_nce)
}

#' Count number of neural network parameters
#' @param weights The weights of a neural network, accessible via model$weights
#' for a tensorflow neural network called model
#' @return num_params (integer): total number of parameters in the network
#' @export
count_NN_params <- function(weights) {
  num_params <- sum(sapply(weights, function(w) {
    shape <- w$shape$as_list()
    prod(unlist(shape))
  }))
  return(num_params)
}

#' Monte Carlo integral of nce_cop (using importance sampling)
#' @param model The classifier trained to distinguish non-simplified data
#' with labels 1 from simplified data with labels 0
#' @param fitted_vine The vine copula fitted to the observed data.
#' @param nu Ratio of noise to true samples
#' @param n_samples Number of samples to create to compute the integral
#' @param data_dim_if_unif Data dimensionality. If this is given (and !=-1),
#' then the samples are drawn from a uniform distribution on [0,1]^d.
#' Otherwise the samples are drawn from the simplified vine
#' @param user_info Whether to display an information, which steps are currently running.
#' @return The Monte Carlo integral approximation.
#' @export
compute_integral <- function(model,
                             fitted_vine,
                             n_samples,
                             nu,
                             data_dim_if_unif=-1,
                             user_info=FALSE){
  samples <- 0
  if (data_dim_if_unif != -1){
    if (user_info){
      print("Start noise sampling (uniform)")
    }
    # simulate random uniform samples on [0,1]^d
    samples <- matrix(runif(n_samples * data_dim_if_unif),
                      ncol=data_dim_if_unif)
  } else {
    if (user_info){
      print("Start noise sampling (simplified vine)")
    }
    # simulate samples from the simplified vine copula
    samples <- rvinecop(n_samples, fitted_vine)
  }
  if(user_info){
    print("Noise samples created")
  }
  p_simp <- 1
  if(data_dim_if_unif==-1){
    # If samples are drawn from the simplified vine copula,
    # the density needs to be evaluated
    if(user_info){
      print("Evaluating noise density (simplified vine)")
    }
    p_simp <- dvinecop(samples, fitted_vine)
  }
  if(user_info){
    print("Evaluating neural network output")
  }
  verbose <- 0
  if(user_info){
    verbose <- 1
  }
  p_nce <- NCE_cop(fitted_vine=fitted_vine, model=model, obs=samples, nu=nu, verbose=verbose)
  return(mean(p_nce/p_simp))
}

#' Compute the g values as defined in the thesis
#' @param model Neural Network for binary classification
#' @param orig_data The observed data.
#' @param nu Defaults to 1. Fraction of noise samples to real samples that
#' was used during training of the neural network.
#' @return A vector containing the log-likelihood ratios between the NCE model and the
#' simplified vine
#' @export
compute_gvals <- function(
    model,
    orig_data,
    nu=1){
  predictions <- model %>% predict(orig_data)
  # here a classifier on the logit level h(u) was trained on the data. to get the values g(u),
  # as defined in the thesis, g(u) = log(nu * p(u)/(1-p(u))) = log(nu * exp(h(u))) = h(u) + log(nu)
  # needs to be computed, where nu is the fraction of number of noise samples to observed samples
  g_vals <- as.numeric(predictions) + log(nu)
  return(g_vals)
}


#' performs a D-Vine quantile regression
#' @param covariates_train matrix or dataframe. The covariates to train the quantile regression method.
#' @param covariates_test matrix or dataframe. The covariates to test the quantile regression method.
#' @param response_train Vector. Values to predict in the quantile regression (train set)
#' @param response_test Vector. Values to predict in the quantile regression (test set)
#' @param family_set_name Defaults to "nonparametric".
#' Argument family_set in the D-vine Quantile regression function.
#' Possible choices include "onepar" and "parametric".
#' @param bottom_quantile_levels Defaults to seq(0.05,0.25,0.05). Vector of lower quantiles,
#' for which it is checked, whether the conditioned quantile estimates are >0.
#' Tests whether the alternative model is better
#' @param top_quantile_levels Defaults to seq(0.75,0.95,0.05). Vector of upper quantiles,
#' for which it is checked, whether the conditioned quantile estimates are <0.
#' Tests whether the simplified model is better
#' @return: A List of 2 vectors, the first containing the number of samples, for which
#' q(bottom_quantile_levels) > 0 holds, the second containing the number of samples, for
#' which q(top_quantile_levels) < 0 holds.
#' @export
perform_quant_reg <- function(
    covariates_train,
    covariates_test,
    response_train,
    response_test,
    family_set_name = "nonparametric",
    bottom_quantile_levels = seq(0.05,0.25,0.05),
    top_quantile_levels = seq(0.75,0.95,0.05)
    ){
  orig_data <- rbind(as.data.frame(covariates_train), as.data.frame(covariates_test))
  cov_train_data <- as.data.frame(covariates_train)
  qreg_data <- cbind(response_train, cov_train_data)
  q_reg <- vinereg::vinereg(response_train ~ . ,
                            family_set=family_set_name,
                            data=qreg_data)
  all_quantile_levels <- c(bottom_quantile_levels, top_quantile_levels)
  train_loss <- rep(0,length(all_quantile_levels))
  test_loss <- rep(0,length(all_quantile_levels))
  pred_train <-predict(q_reg, newdata=covariates_train, alpha=all_quantile_levels)
  pred_test <- predict(q_reg, newdata=covariates_test, alpha=all_quantile_levels)
  for(i in 1:length(all_quantile_levels)){
    train_loss[i] <- pinball_loss(y_true=response_train,
                                  y_pred=pred_train[,i],
                                  tau=all_quantile_levels[i])
    test_loss[i] <- pinball_loss(y_true=response_test,
                                 y_pred=pred_test[,i],
                                 tau=all_quantile_levels[i])
  }
  bottom_quantiles <- predict(q_reg, newdata=orig_data, alpha=bottom_quantile_levels)
  top_quantiles <- predict(q_reg, newdata=orig_data, alpha=top_quantile_levels)
  alternative_better <- rep(0,length(bottom_quantile_levels))
  simp_better <- rep(0,length(top_quantile_levels))
  for(i in 1:nrow(orig_data)){
    for(j in 1:length(bottom_quantile_levels)){
      if(bottom_quantiles[i,j] > 0){
        alternative_better[j] <- alternative_better[j] + 1
      }
    }
    for(j in 1:length(top_quantile_levels)){
        if(top_quantiles[i,j] < 0){
          simp_better[j] <- simp_better[j] + 1
        }
    }
  }
  output <- list(alternative_better,
                 simp_better,
                 train_loss,
                 test_loss
                 )
  return(output)
}

#' performs a linear quantile regression
#' @param covariates_train matrix or dataframe. The covariates to train the quantile regression method.
#' @param covariates_test matrix or dataframe. The covariates to test the quantile regression method.
#' @param response_train Vector. Values to predict in the quantile regression (train set)
#' @param response_test Vector. Values to predict in the quantile regression (test set)
#' @param bottom_quantiles_lin Defaults to c(0.05,0.1). Vector of lower quantiles,
#' for which it is checked, whether the conditioned quantile estimates are >0.
#' Tests whether the alternative model is better
#' @param top_quantiles_lin Defaults to c(0.05,0.1). Vector of upper quantiles,
#' for which it is checked, whether the conditioned quantile estimates are <0.
#' Tests whether the simplified model is better
#' @param method Defaults to "br". Which method to use for the linear quantile regression
#' @return: A List of 4 vectors, the first containing the number of samples, for which
#' q(bottom_quantile_levels) > 0 holds, the second containing the number of samples, for
#' which q(top_quantile_levels) < 0 holds. Behind that train and test pinball loss are listed.
#' @export
perform_linear_quant_reg <- function(
    covariates_train,
    covariates_test,
    response_train,
    response_test,
    bottom_quantiles_lin=seq(0.05,0.25,0.05),
    top_quantiles_lin=seq(0.75,0.95,0.05),
    method="br"){
  orig_data <- rbind(as.data.frame(covariates_train), as.data.frame(covariates_test))
  # perform the quantile regression
  linear_qreg_dataframe <- as.data.frame(covariates_train)
  all_quantiles <- c(bottom_quantiles_lin, top_quantiles_lin)
  fit_linear_model <- quantreg::rq(response_train ~ .,
                         data=linear_qreg_dataframe,
                         tau=all_quantiles)
  # calculate training and test pinball loss
  train_loss <- rep(0,length(all_quantiles))
  test_loss <- rep(0, length(all_quantiles))
  lin_pred_train <-predict(fit_linear_model,
                           newdata=as.data.frame(covariates_train))
  lin_pred_test <- predict(fit_linear_model,
                           newdata=as.data.frame(covariates_test))
  for(i in 1:length(all_quantiles)){
    train_loss[i] <- pinball_loss(y_true=response_train,
                                  y_pred=lin_pred_train[,i],
                                  tau=all_quantiles[i])
    test_loss[i] <- pinball_loss(y_true=response_test,
                                 y_pred=lin_pred_test[,i],
                                 tau=all_quantiles[i])
  }
  # create the output of the relevant fields
  preds_lin <- predict(fit_linear_model,
                       newdata = as.data.frame(orig_data))
  alternative_better_lin <- rep(0, length(bottom_quantiles_lin))
  simp_better_lin <- rep(0, length(top_quantiles_lin))
  for(i in 1:nrow(orig_data)){
    for(j in 1:length(bottom_quantiles_lin)){
      if(preds_lin[i,j] > 0){
        alternative_better_lin[j] <- alternative_better_lin[j] + 1
      }
    }
    for(j in 1:length(top_quantiles_lin)){
      if(preds_lin[i,(j+length(bottom_quantiles_lin))] < 0){
        simp_better_lin[j] <- simp_better_lin[j] + 1
      }
    }
  }
  return(list(alternative_better_lin, simp_better_lin, train_loss, test_loss))
}

#' performs a Neural Network quantile regression (MCQRNN)
#' @param covariates_train matrix or dataframe. The covariates to train the quantile regression method.
#' @param covariates_test matrix or dataframe. The covariates to test the quantile regression method.
#' @param response_train Vector. Values to predict in the quantile regression (train set)
#' @param response_test Vector. Values to predict in the quantile regression (test set)
#' @param bottom_q_NN_train Vector of lower quantiles, for which the model is fitted
#' Note: a longer vector here tends to improve the fit but also increase runtime considerably.
#' Recommendation: Train on no more than every 5 percent quantiles, predict the others with mcqrnn
#' @param bottom_q_NN_predict Vector of lower quantiles,
#' for which it is checked, whether the conditioned quantile estimates are >0.
#' Tests whether the alternative model is better.
#' @param top_q_NN_train Vector of top quantiles, for which the model is fitted
#' Note: a longer vector here tends to improve the fit but also increase runtime considerably.
#' Recommendation: Train on no more than every 5 percent quantiles, predict the others with mcqrnn
#' @param top_q_NN_predict Vector of top quantiles,
#' for which it is checked, whether the conditioned quantile estimates are >0.
#' Tests whether the simplified model is better.
#' @param num_hidden int, defaults to 10. Number of hidden units in the
#' single hidden layer of the NN
#' @param num_trials int, defaults to 1. Number of repeated fitting procedures
#' (To try and avoid local minima)
#' @param penalty float, defaults to 0.1. Parameter for weight decay regularization
#' @param max_iter int, defaults to 500. Maximum number of iterations of the
#' optimization algorithm
#' @param activation Which activation function to use. Needs to be one of the activation
#' functions in the qrnn package. Defaults to sigmoid.
#' @param user_info Boolean, defaults to FALSE. If True, diagnostic messages are
#' printed during optimization
#' @return A List of 4 vectors. The first vector contains the number of samples, for which
#' q(bottom_quantile_levels) > 0 holds, the second containing the number of samples, for
#' which q(top_quantile_levels) < 0 holds. Behind that train and test pinball loss are listed.
#' @export
perform_quant_reg_mcqrnn <- function(
    covariates_train,
    covariates_test,
    response_train,
    response_test,
    bottom_q_NN_train = c(0.01,0.05,0.1,0.15,0.2),
    bottom_q_NN_predict=seq(0.01,0.2,0.01),
    top_q_NN_train=c(0.8,0.85,0.9,0.95,0.99),
    top_q_NN_predict = seq(0.8,0.99,0.01),
    num_hidden = 10,
    num_trials=1,
    penalty=0.1,
    max_iter=500,
    activation = qrnn::sigmoid,
    user_info=FALSE
    ){
  orig_data <- as.matrix(rbind(covariates_train, covariates_test))
  # perform quantile regression
  q_NN_train <- c(bottom_q_NN_train, top_q_NN_train)
  x_qrnn <- as.matrix(covariates_train)
  y_qrnn <- matrix(response_train,ncol=1)
  fitted_mcqrnn <- qrnn::mcqrnn.fit(x=x_qrnn,
                              y=y_qrnn,
                              tau=q_NN_train,
                              n.hidden=num_hidden,
                              n.trials=num_trials,
                              penalty=penalty,
                              Th=activation,
                              iter.max=max_iter,
                              trace=user_info)
  # calculate training and test pinball loss
  train_loss <- rep(0,length(q_NN_train))
  test_loss <- rep(0, length(q_NN_train))
  q_NN_predict <- c(bottom_q_NN_predict, top_q_NN_predict)
  qrnn_pred_train <-qrnn::mcqrnn.predict(covariates_train, fitted_mcqrnn, tau=q_NN_predict)
  qrnn_pred_test <- qrnn::mcqrnn.predict(covariates_test, fitted_mcqrnn, tau=q_NN_predict)
  for(i in 1:length(q_NN_predict)){
    train_loss[i] <- pinball_loss(y_true=response_train,
                               y_pred=qrnn_pred_train[,i],
                               tau=q_NN_predict[i])
    test_loss[i] <- pinball_loss(y_true=response_test,
                              y_pred=qrnn_pred_test[,i],
                              tau=q_NN_predict[i])
  }
  qrnn_preds <- qrnn::mcqrnn.predict(orig_data, fitted_mcqrnn, tau=q_NN_predict)
  alternative_better_qrnn <- rep(0, length(bottom_q_NN_predict))
  simp_better_qrnn <- rep(0, length(top_q_NN_predict))
  for(i in 1:nrow(orig_data)){
    for(j in 1:length(bottom_q_NN_predict)){
      if(qrnn_preds[i,j] > 0){
        alternative_better_qrnn[j] <- alternative_better_qrnn[j] + 1
      }
    }
    for(j in 1:length(top_q_NN_predict)){
      if(qrnn_preds[i,(j+length(bottom_q_NN_predict))] < 0){
        simp_better_qrnn[j] <- simp_better_qrnn[j] + 1
      }
    }
  }
  return(list(alternative_better_qrnn, simp_better_qrnn, train_loss, test_loss))
}

#' Calculate the pinball loss (sometimes called check loss)
#' @param y_true observed values
#' @param y_pred predicted quantiles
#' @param tau quantile level (in the thesis this is alpha,
#' tau for consistency with the other packages)
#' @return A number, which is the pinball loss for the given quantile level
#' alpha, true values y_true and predicted values y_pred
#' @export
pinball_loss <- function(y_true, y_pred, tau) {
  differences <- y_true - y_pred
  return(mean(pmax(tau * differences, (tau - 1) * differences)))
}

#' Train Test split for the data used in quantile regression
#' @param covariates matrix. The covariates, used to derive conditional quantiles
#' @param predictors vector or matrix. The values for which quantiles should be estimated
#' @param train_proportion Number, defaults to 0.9. Fraction of the data to use for training
#' @return list with 4 elements: covariates_train, covariates_test, predictors_train, predictors_test
#' (i.e. original covariates and predictors split up randomly)
#' @export
qreg_train_test_split <- function(covariates, predictors, train_proportion){
  num_train <- floor(nrow(covariates) * train_proportion)
  train_indices <- sample(1:nrow(covariates), size = num_train)
  covariates_train <- covariates[train_indices,]
  covariates_test <- covariates[-train_indices,]
  predictors_train <- predictors[train_indices]
  predictors_test <- predictors[-train_indices]
  return(list(covariates_train, covariates_test, predictors_train, predictors_test))
}

# Example Usage:
# split_output <- train_test_split(orig_data, simplified_data)
# x_train <- split_output[[1]]
# x_test <- split_output[[2]]
# y_train <- split_output[[3]]
# y_test <- split_output[[4]]
# model <- build_model(input_dim=ncol(x_train))
# model %>% summary
# model <- train_model(model, x_train, y_train, num_epochs=10)
# model %>% evaluate(x_test, y_test)
# g_vals <- compute_gvals(model, orig_data)
# output <- perform_quant_reg(g_vals, orig_data)
# print(output)
