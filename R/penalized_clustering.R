
# Function for fitting penalized MBC models \citep{Zhou2009}, with groupwise penalizations---------------------
# the penalization type is set via the group_shrinkage function argument:
# - "common": Zhou2009 model
# - "weighted_by_W0": 1/omega_k^(0) used as P_k
# - "weighted_by_dist_to_I": distance bw omega_k^(0) and diagonal matrix used as P_k (all elements of P_k are equal in this case)
# - "weighted_by_cellwise_dist_to_I": cellwise distance bw omega_k^(0) and diagonal matrix used as P_k
#' @export
fit_penalized_clust <-
  function(data,
           K,
           lambda_omega, # common lambda for all groups (it acts as a multiplicative factor for P_k)
           lambda_mu,
           penalize_diag=FALSE,
           group_shrinkage = c(
             "common",
             "weighted_by_W0",
             "weighted_by_dist_to_I",
             "weighted_by_cellwise_dist_to_I"
           ),
           initialization=NULL,
           control_EM_algorithm =control_EM()
  ) {


    data <- data.matrix(scale(data)) # work with standardized data

    N <- nrow(data)
    p <- ncol(data)

    #### EM ######################################################################


    if(length(lambda_mu)==1){
      lambda_mu <- rep(lambda_mu, p)
    }
    if (!is.matrix(lambda_omega))
      lambda_omega <- matrix(lambda_omega, nrow = p, ncol = p)
    if (penalize_diag == FALSE) diag(lambda_omega) <- 0

    # intialization of cluster allocation
    if(is.null(initialization)){
      hc_init <- mclust::hc(data, modelName = "VVV", use = "SVD")
      z <- mclust::unmap(mclust::hclass(hc_init, K))
    } else {
      z <- mclust::unmap(initialization)
    }

    temp <- mclust::covw(data, z, normalize = FALSE)

    ######## I set the group-wise penalization via P_k

    P_k <- switch (
      group_shrinkage,
      "common" = array(1, dim = c(p, p, K)),
      "weighted_by_W0" = 1 / abs(temp$W),
      "weighted_by_dist_to_I" = 1 / array(data = rep(
        apply(temp$W, 3, function(omega)
          shapes::distcov(S1 = omega, S2 = diag(p), method = "Riemannian")), each =p ^ 2),
        dim = c(p, p, K)), # I use Riemannian metric, but since PD space is non-euclidian any distance in Dryden2009 can be employed
      "weighted_by_cellwise_dist_to_I" = 1/array(apply(temp$W, 3, function(omega)
        abs(omega - diag(
          p
        ))), dim = c(p, p, K))
    )

    # start EM parameters and containers

    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err

    iter <- 0
    loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
    loglik_pen_vec <- NULL
    penalty <- penalty_mu <- rep(0, K)
    mu <- matrix(NA, p, K)
    sigma <- omega <- array(NA, c(p, p, K))

    for (k in 1:K) {
      omega[, , k] <- temp$W[, , k] # starting value for omega
    }

    # dimnames(sigma) <- dimnames(omega) <- list(varnames, varnames)
    crit <- TRUE


    # run EM
    while (crit) {
      #### M step -------------------------------------------------
      nk <- colSums(z)
      pro <- nk / N

      if(iter!=0){ # I do not do this the first iteration as I have already computed it outside the while loop
        temp <- mclust::covw(data, z)  # compute weighted quantities
      }


      mu <- temp$mean    # sample mean

      # Sparse mu estimation via coordinate ascent algorithm
      if (any(lambda_mu != 0)) {
        for (k in 1:K) {
          for (j in 1:p) {
            first_piece_check <-
              c(sweep(
                x = data[, -j],
                MARGIN = 2,
                STATS = mu[-j, k],
                FUN = "-"
              ) %*% omega[-j, j, k])

            second_piece_check <- data[, j] * omega[j, j, k]

            check_sparse_mu <-
              abs(sum(z[, k] * (
                first_piece_check + second_piece_check
              ))) < c(lambda_mu)[j]

            if (check_sparse_mu) {
              mu[j, k] <- 0
            } else {
              numerator_penalized_mu <-
                sum(z[, k] * c(data %*% omega[, j, k])) - nk[k] *
                (c(mu[, k] %*% omega[, j, k]) - mu[j, k] * omega[j, j, k])

              mu[j, k] <-
                (numerator_penalized_mu - sign(numerator_penalized_mu) * c(lambda_mu)[j]) /
                (nk[k] * omega[j, j, k])
            }
          }
          penalty_mu[k] <- sum(abs(lambda_mu * mu[, k]))
        }
      }

      # set initial covariance/precision matrix for glasso - ensure monotonicity
      if (iter < 1) {
        start <- "cold"
        start_sigma <- temp$S
        start_omega <- array(apply(temp$S, 3, solve), c(p, p, K))
      } else {
        start <- "warm"
        start_sigma <- sigma
        start_omega <- omega
      }

      if (any(lambda_omega != 0)) {
        for (k in 1:K) {
          # graphical lasso estimation
          gl <- glasso::glasso(
            temp$S[, , k],
            rho = 2 * lambda_omega * P_k[, , k] / nk[k],
            start = start,
            w.init = start_sigma[, , k],
            wi.init = start_omega[, , k]
          )
          sigma[, , k] <- gl$w
          omega[, , k] <- gl$wi
          penalty[k] <-
            sum(abs(lambda_omega * P_k[, , k] * omega[, , k]))
        }
      } else {
        sigma <- temp$S
        omega <- temp$W
      }

      ####---------------------------------------------------------

      #### E step -------------------------------------------------
      dens <- matrix(NA, N, K)
      for (k in 1:K)
        dens[, k] <- mclust::dmvnorm(data, mu[, k], sigma[, , k], log = TRUE)
      denspro <- sweep(dens, 2, log(pro), "+")

      # zetas
      zMax <- apply(denspro, 1, max)
      loghood <- zMax + log(rowSums(exp(denspro - zMax)))
      z <- exp(denspro - loghood)
      ####---------------------------------------------------------

      #### loglik -------------------------------------------------
      loglik <- sum(loghood)
      loglik_pen <- loglik - sum(penalty) - sum(penalty_mu)

      # check convergence
      err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
      loglik_pen_prev <- loglik_pen
      loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
      iter <- iter + 1

      crit <- (err > tol & iter < itermax)
      #------------------------------------------------------------
    }

    # Results collection

    parameters <- list(pro=pro, mu=mu, sigma=sigma, omega=omega)
    n_par_pro <- K - 1
    n_par_mu <- sum(!mu == 0)
    n_par_omega <- p * K + sum(apply(omega,
                                     3, function(A) A[upper.tri(A)] != 0))
    bic_final <- 2 * loglik - (n_par_pro + n_par_mu + n_par_omega) * log(N)

    OUT <-
      list(
        loglik = loglik,
        loglik_pen = loglik_pen,
        parameters = parameters,
        z = z,
        classification = mclust::map(z),
        bic = bic_final,
        n_params = c(pro = n_par_pro, mu = n_par_mu, omega = n_par_omega),
        penalty = list(lambda_mu = lambda_mu,
                       lambda_omega = lambda_omega,
                       P_k = P_k),
        LLK_trace = loglik_pen_vec,
        iter = iter
      )
    OUT
  }
