# generate_data_vanilla <- function (n, # sample size
#                                    a, # intercept
#                                    b,
#                                    gamma, # deterministic coefficient
#                                    p,
#                                    dgp_para = 2) {
#
#     chisq_df <- dgp_para
#
#     group_id <- as.integer(runif(n) > p) + 1
#     b0 <- b[group_id]
#
#     x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
#     z_1 <- x + rnorm(n)
#     z_2 <- z_1 + rnorm(n)
#     z <- cbind(z_1, z_2)
#     sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
#     u <- rnorm(n, mean = 0, sd = sqrt(sigma_i))
#
#     y <- a + x * b0 + as.numeric(z %*% gamma) + u
#
#     return(cbind(y, x, z))
# }
#
# generate_data_vanilla_3 <- function(n, # sample size
#                                     a, # intercept
#                                     b,
#                                     gamma, # deterministic coefficient
#                                     p,
#                                     dgp_para = 2) {
#     chisq_df <- dgp_para
#
#     sample_ind <- runif(n)
#     group_id <- rowSums(
#         sapply(p, function(x){
#             as.integer(sample_ind > x)
#         })
#     )
#     b0 <- b[group_id]
#
#     x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
#     z_1 <- x + rnorm(n)
#     z_2 <- z_1 + rnorm(n)
#     z <- cbind(z_1, z_2)
#     sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
#     u <- rnorm(n, mean = 0, sd = sqrt(sigma_i))
#
#     y <- a + x * b0 + as.numeric(z %*% gamma) + u
#
#     return(cbind(y, x, z))
# }

# set.seed(10000)
# n <- 1000
# a <- 0.5
# b <- c(1, 2, 3)
# p <- c(0, 0.25, 0.75)
# gamma0 <- c(1, 1)
#
# for (r in 1:1000) {
#     t0 <- Sys.time()
#     d <- generate_data_vanilla_3(n, a, b, gamma0, p)
#     x <- d[, 'x']
#     y <- d[, 'y']
#     z <- d[, -c(1, 2)]
#
#     res3 <- ccrm_est_K3(x, y, z, s_max = 6)
#     print(r)
#     print(Sys.time() - t0)
# }
