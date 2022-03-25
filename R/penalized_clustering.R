# WRAPPER FUNCTION
# Function for fitting penalized MBC models \citep{Zhou2009}, with groupwise penalizations---------------------
# the penalization type is set via the group_shrinkage function argument:
# - "common": Zhou2009 model
# - "weighted_by_W0": 1/omega_k^(0) used as P_k
# - "weighted_by_dist_to_I": distance bw omega_k^(0) and diagonal matrix used as P_k (all elements of P_k are equal in this case)
# - "weighted_by_dist_to_diag_W0": distance bw omega_k^(0) and diagonal W0 used as P_k (all elements of P_k are equal in this case)
# - for the last two group_shrinkage type, given that the PD space is non-Euclidian, several distances may be considered, controlled with the argument
# distance_method=c("Procrustes","ProcrustesShape","Riemannian","Cholesky", "Euclidean", "LogEuclidean", "RiemannianLe"). See Dryden 2009 AOAS for a full account on the problem.
#' @export
fit_penalized_clust <-
  function(data,
           K,
           lambda_omega,
           # common lambda for all groups (it acts as a multiplicative factor for P_k)
           lambda_mu,
           penalize_diag = FALSE,
           group_shrinkage = c(
             "common",
             "weighted_by_W0",
             "weighted_by_dist_to_I",
             "weighted_by_dist_to_diag_W0",
           ),
           distance_method = "Euclidean",
           # used in "weighted_by_dist_to_I" and "weighted_by_dist_to_diag_W0" only
           initialization = NULL,
           lambda_omega_0 = 50,
           # when the dimensionality is high, Omega_0 might be singular, so we use glasso with penalty equal lambda_omega_0 to start the procedure
           epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
           # the denominator of formula (6)
           control_EM_algorithm = control_EM(),
           verbose = interactive()) {
    # Wrapper function for performing penalized MBC varying hyper-parameters
    # The best model is the one that maximizes the BIC
    start_time <- Sys.time()
    # data <- data.matrix(scale(data))
    data <- data.matrix(data)

    N <- nrow(data)
    p <- ncol(data)

    all_hyperparameters <-
      expand.grid(K = K,
                  lambda_omega = lambda_omega,
                  lambda_mu = lambda_mu)
    n_different_models <- nrow(all_hyperparameters)

    models_container <-
      vector(mode = "list", length = n_different_models)

    if (verbose) {
      cat("Fitting:\n")
      utils::flush.console()
      pbar <- utils::txtProgressBar(min = 0,
                                    max = n_different_models,
                                    style = 3)
      on.exit(close(pbar))
      ipbar <- 0
    }

    # Initialization common for all models
    if (is.null(initialization)) {
      global_init <- mclust::hc(data, modelName = "VVV", use = "SVD")
      # z <- mclust::unmap(mclust::hclass(hc_init, K)) # FIXME to be removed
    } else {
      global_init <- initialization
    }

    for (model in 1:n_different_models) {
      models_container[[model]] <-
        fit_single_parameters_set(
          data = data,
          K = all_hyperparameters[model, "K"],
          lambda_omega = all_hyperparameters[model, "lambda_omega"],
          lambda_mu = all_hyperparameters[model, "lambda_mu"],
          penalize_diag = penalize_diag,
          group_shrinkage = group_shrinkage,
          distance_method = distance_method,
          initialization = global_init,
          lambda_omega_0 = lambda_omega_0,
          epsilon_weighted_by_W0 = epsilon_weighted_by_W0,
          control_EM_algorithm = control_EM_algorithm,
          N = N,
          p = p
        )

      if (verbose) {
        ipbar <- ipbar + 1
        utils::setTxtProgressBar(pbar, ipbar)
      }
    }

    models_BICS <- sapply(models_container, "[[", "bic")
    max_bic_model <- which.max(models_BICS)

    selected_model <- models_container[[max_bic_model]]
    selected_model$BIC <- cbind(all_hyperparameters, bic=models_BICS)

    end_time <- Sys.time()
    selected_model$elapsed_time <- end_time - start_time
    if (verbose) {
      cat(
        "\nModel with K=",
        selected_model$K,
        ", lambda_omega=",
        all_hyperparameters[max_bic_model, "lambda_omega"],
        ", lambda_mu=",
        all_hyperparameters[max_bic_model, "lambda_mu"],
        " returns the highest BIC.",
        sep = ""
      )
    }
    selected_model
  }
