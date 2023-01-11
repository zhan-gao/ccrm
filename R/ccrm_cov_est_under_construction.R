init_moment_est <- function(x, y, z, K) {

  require("multicool")
  z <- cbind(1, z)
  R <- 2 * K - 1

  # store Eb_0 and Eu_0 to be 1 as the first element
  Eb_vec <- c(1, rep(0, R))
  Eu_vec <- c(1, rep(0, R))

  # OLS estimation
  w <- cbind(x, z)
  phi_hat <- solve(t(w) %*% w) %*% (t(w) %*% y)

  gamma_hat <- phi_hat[-1]
  Eb_vec[2] <- phi_hat[1] # first element to be 1

  for (r in 2:R) {

    if (r %% 2 == 0) {
      s = 2
    } else {
      s = 1
    }

    # Get all combinations of q_1, q_2, q_3
    if (r == 2) {
      multinom_order <- lapply(list(c(2,0,0), c(1,1,0)), function(l){allPerm(initMC(l))})
    } else {
      multinom_order <- lapply(genComp(r, 3, TRUE), function(l){allPerm(initMC(l))})
    }
    multinom_order <- do.call(rbind, multinom_order)
    # exclude EB_r and Eu_r
    multinom_order <- multinom_order[(multinom_order[, 3] != r) & (multinom_order[, 1] != r), ]

    # Compute multinomial coefficients
    multinom_coef <- apply(multinom_order, 1, function(l){multinom(l, counts = TRUE)})

    # compute the summation term (by a loop)
    num_term <- length(multinom_coef)
    sum_term_vec <- rep(0, num_term)
    sum_term_s_vec <- rep(0, num_term)

    for (q in 1:num_term) {

      order_q <- multinom_order[q, ]
      coef_q <- multinom_coef[q]

      E_xz <- mean(x ^ order_q[1] * (z %*% gamma_hat) ^ order_q[2])
      E_xz_s <- mean(x ^ (order_q[1] + s) * (z %*% gamma_hat) ^ order_q[2])

      sum_term_vec[q] <- coef_q * E_xz * Eb_vec[order_q[1] + 1] * Eu_vec[order_q[3] + 1]
      sum_term_s_vec[q] <- coef_q * E_xz_s * Eb_vec[order_q[1] + 1] * Eu_vec[order_q[3] + 1]
    }

    # Compute Eb_r
    Eb_vec[r + 1] <- ((mean(x^s * y^r) - mean(x^s) * mean(y^r)) -
                  (sum(sum_term_s_vec) - sum(sum_term_vec) * mean(x^s))) /
                 (mean(x^(r + s)) - mean(x^r) * mean(x^s))

    # Compute the error moment
    Eu_vec[r + 1] <- mean(y ^ r) - sum(sum_term_vec) - mean(x ^ r) * Eb_vec[r + 1]
  }

  list(
    gamma = gamma_hat,
    Eb_hat = Eb_vec[-1]
  )
}



multinom_gen <- function(r, exclude = TRUE) {

  if (r <= 1) {
    return(list(order = t(rep(0,3)),
                coef = 1))
  }

  # Get all combinations of q_1, q_2, q_3
  if (r == 2) {
    multinom_order <- lapply(list(c(2,0,0), c(1,1,0)), function(l){allPerm(initMC(l))})
  } else {
    multinom_order <- lapply(genComp(r, 3, TRUE), function(l){allPerm(initMC(l))})
  }
  multinom_order <- do.call(rbind, multinom_order)

  if (exclude) {
    multinom_order <- multinom_order[(multinom_order[, 3] != 1), ]
  }
  # Compute multinomial coefficients
  multinom_coef <- apply(multinom_order, 1, function(l){multinom(l, counts = TRUE)})


  list(order = multinom_order,
       coef = multinom_coef)
}

gmm_moment_est <- function(x, y, z, K, s_max = 2 * K) {

  require(multicool)

  n <- nrow(z)
  R <- 2 * K - 1

  # initial moment estimator
  init_est_result <- init_moment_est(x, y, z, K)
  theta_init_est <- c(init_est_result$gamma[1], init_est_result$Eb_hat)
  # Get rid of z * gamma
  y <- y - as.numeric(z %*% init_est_result$gamma[-1])

  # ----- Moment conditions



}


error_moment <- function(x, y, z, gamma, Eb, r_max) {

  # don't want this.

  # r_max >= 3
  require("multicool")

  Eu = rep(0, r_max) # Note that we assume u is mean 0

  for(r in 2:r_max){

    # Get all combinations of q_1, q_2, q_3
    if (r == 2) {
      multinom_order <- lapply(list(c(2,0,0), c(1,1,0)), function(l){allPerm(initMC(l))})
    } else {
      multinom_order <- lapply(genComp(r, 3, TRUE), function(l){allPerm(initMC(l))})
    }
    multinom_order <- do.call(rbind, multinom_order)
    multinom_order <- multinom_order[multinom_order[, 3] != r, ]

    # Compute multinomial coefficients
    multinom_coef <- apply(multinom_order, 1, function(l){multinom(l, counts = TRUE)})

    # compute the summation term (by a loop)
    num_term <- length(multinom_coef)
    sum_term_vec <- rep(0, num_term)


    for (q in 1:num_term) {
      order_q <- multinom_order[q, ]
      coef_q <- multinom_coef[q]
      E_xz <- mean(x^order_q[1] * (z%*%gamma)^order_q[2])
      sum_term_vec[q] <- coef_q * E_xz * Eb[order_q[1]] * Eu(order_q[3])
    }

    # Compute the error moment
    Eu[r] <- mean(y^r) - sum(sum_term_vec)
  }

  Eu
}



