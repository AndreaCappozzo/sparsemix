
fit_single_parameters_set <-
  function(data,
           K,
           lambda_omega,
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
           N,
           p) {


    #### EM ######################################################################


    if(length(lambda_mu)==1) {
      lambda_mu <- rep(lambda_mu, p)
    }

    if (!is.matrix(lambda_omega)) {
      lambda_omega <- matrix(lambda_omega, nrow = p, ncol = p)
    }

    if (!is.matrix(lambda_omega_0)) {
      lambda_omega_0 <- matrix(lambda_omega_0, nrow = p, ncol = p)
    }

    if (penalize_diag == FALSE) {
      diag(lambda_omega) <- 0
      diag(lambda_omega_0) <- 0
    }

    # initialization of cluster allocation

    if(class(initialization)=="hc"){
      z <- mclust::unmap(mclust::hclass(initialization, K))
    } else {
      z <- mclust::unmap(initialization)
    }

    temp <- mclust::covw(data, z, normalize = FALSE)

    # Check if I have issues in the initial estimation of S
    # determinant_S <- apply(temp$S, 3, function(X)
    #   determinant(X))
    # log_det_S <- sapply(determinant_S, FUN = "[[", 1)
    # sign_det_S <- sapply(determinant_S, FUN = "[[", 2)
    #
    # singularity_problem <-
    #   (any(sign_det_S<0))|any(log_det_S<log(.Machine$double.neg.eps)) # I check if the matrices temp$S are invertible
    #
    # alternative
    eig <- sapply(1:K, function(k) eigen(temp$S[,,k], only.values = TRUE)$val)
    singularity_problem <- any(eig < .Machine$double.eps)

    if (!singularity_problem) {
      # if temp$S are invertible, then I compute the Omega_0 as the inverse of S_0
      omega_0 <- array(apply(temp$S, 3, solve), dim = c(p, p, K))
      sigma_0 <- temp$S
    } else{
      # if temp$S are singular, I use graphical lasso to compute the initial estimates
      omega_0 <- sigma_0 <- array(dim = c(p, p, K))
      nk <- colSums(z)
      for (k in 1:K) {
        # graphical lasso estimation
        gl <- glassoFast::glassoFast(S = temp$S[, , k],
                                     rho = 2 * lambda_omega_0 / nk[k])
        omega_0[, , k] <- gl$wi
        sigma_0[, , k] <- gl$w
      }
    }
    ######## I set the group-wise penalization via P_k

    P_k <- switch (
      group_shrinkage,
      "common" = array(1, dim = c(p, p, K)),
      "weighted_by_W0" = 1 / (epsilon_weighted_by_W0 + abs(omega_0)),
      "weighted_by_dist_to_I" = 1 / array(data = rep(
        apply(omega_0, 3, function(omega)
          shapes::distcov(S1 = omega, S2 = diag(p), method = distance_method)), each =p ^ 2),
        dim = c(p, p, K)),
      "weighted_by_dist_to_diag_W0" = 1 / array(data = rep(
        apply(omega_0, 3, function(omega)
          shapes::distcov(S1 = omega, S2 = diag(diag(omega)), method = distance_method)), each =p ^ 2),
        dim = c(p, p, K))
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

    # start omega and sigma
    omega <- omega_0
    sigma <- sigma_0


    # dimnames(sigma) <- dimnames(omega) <- list(varnames, varnames)
    crit <- TRUE
    z_prev <- z
    to_break <- FALSE

    # ME algorithm

    while (crit) {

      #### M step -------------------------------------------------
      nk <- colSums(z)
      pro <- nk / N

      if(iter != 0){ # I do not do this the first iteration as I have already computed it outside the while loop
        temp <- mclust::covw(data, z, normalize=FALSE)  # compute weighted quantities
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

          all_zeros <- apply(mu, 2, function(mu_g) all(mu_g == 0))

          if (any(all_zeros)) {
            warning("An entire mean component vector was shrunk to 0: consider decreasing lambda_mu")
          }
          penalty_mu[k] <- sum(abs(lambda_mu * mu[, k]))
        }
      }

      if ( any(lambda_omega != 0) ) {
        for (k in 1:K) {
          # graphical lasso estimation
          # gl <- glassoFast::glassoFast(
          #   S=temp$S[, , k],
          #   rho = 2 * lambda_omega * P_k[, , k] / nk[k]
          # )

          gl <- try(glassoFast::glassoFast(S = temp$S[,,k],
                                           rho = 2 * lambda_omega * P_k[, , k] / nk[k],
                                           start = "warm",
                                           wi.init = omega[,,k],
                                           w.init = sigma[, , k]), silent = TRUE)
          if ( class(gl) == "try-error" | any(is.nan(gl$wi)) ) {
            gl <- glassoFast::glassoFast(S = temp$S[,,k],
                                         rho = 2 * lambda_omega * P_k[, , k] / nk[k])
          }

          if( any( diag(gl$w) < sqrt(.Machine$double.eps) ) ) {
            to_break <- TRUE
            break
          }

          sigma[, , k] <- gl$w
          omega[, , k] <- gl$wi
          penalty[k] <-
            sum(abs(lambda_omega * P_k[, , k] * omega[, , k]))
        }
      } else {
        sigma <- temp$S
        omega <- array(apply(temp$S, 3, solve), c(p, p, K))
      }


      # check for near zero variance and break
      if ( to_break ) {
        if ( iter != 0 ) {

          nk <- colSums(z_prev)
          pro <- nk / N
          temp <- mclust::covw(data, z_prev, normalize = FALSE)
          mu <- temp$mean
          # TODO: fix mean estimation to account for sparsity

          if (any(lambda_omega != 0)) {

            for (k in 1:K) {
              gl <- try(glassoFast::glassoFast(S = temp$S[,,k],
                                               rho = 2 * lambda_omega * P_k[, , k] / nk[k],
                                               start = "warm",
                                               wi.init = omega[,,k],
                                               w.init = sigma[, , k]), silent = TRUE)
              if ( class(gl) == "try-error" | any(is.nan(gl$wi)) ) {
                gl <- glassoFast::glassoFast(S = temp$S[,,k],
                                             rho = 2 * lambda_omega * P_k[, , k] / nk[k])
              }

              sigma[, , k] <- gl$w
              omega[, , k] <- gl$wi
              penalty[k] <-
                sum(abs(lambda_omega * P_k[, , k] * omega[, , k]))
            }
          } else {
            sigma <- temp$S
            omega <- array(apply(temp$S, 3, solve), c(p, p, K))
          }

        } else {
          OUT <- list(loglik = NA, loglik_pen = NA, parameters = NA, z = z, K = K, classification = NA,
                      bic = NA, n_params = NA, penalty = list(lambda_mu = lambda_mu, lambda_omega = lambda_omega, P_k = P_k),
                      LLK_trace = NA, iter = iter)

          return( OUT )
        }
      }


      #### E step -------------------------------------------------

      dens <- matrix(NA, N, K)
      for (k in 1:K){
        dens[, k] <- mclust::dmvnorm(data, mu[, k], sigma[, , k], log = TRUE)
      }
      denspro <- sweep(dens, 2, log(pro), "+")

      # zetas
      zMax <- apply(denspro, 1, max)
      loghood <- zMax + log(rowSums(exp(denspro - zMax)))
      z_prev <- z
      z <- exp(denspro - loghood)

      #### loglik -------------------------------------------------
      loglik <- sum(loghood)
      loglik_pen <- loglik - sum(penalty) - sum(penalty_mu)

      # check convergence
      err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
      loglik_pen_prev <- loglik_pen
      loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
      iter <- iter + 1

      crit <- if (to_break) FALSE else (err > tol & iter < itermax)
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
        K=K,
        classification = mclust::map(z),
        bic = bic_final,
        n_params = c(pro = n_par_pro, mu = n_par_mu, omega = n_par_omega),
        penalty = list(lambda_mu = lambda_mu,
                       lambda_omega = lambda_omega,
                       P_k = P_k),
        LLK_trace = loglik_pen_vec,
        iter = iter
      )

    return( OUT )
  }
