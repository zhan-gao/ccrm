#' The wrapper of estimation functions
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-1)
#' @param z control varibles (N-by-p)
#' @param theta_init initial value for distribution parameters
#' @param s_max
#' @param est_method
#'
#' @return A list contains estimated coefficients and inferential statistics
#' \item{theta_b}{Estimated distributional parameter p, b_L, b_H}
#' \item{theta_m}{Estimated moments of beta_i}
#' \item{theta_u}{Estimated moments of u_i}
#' \item{gamma}{estimated coefficients of control varibles}
#' \item{V_theta}{variance}
#'
#' @export
#'
#'
ccrm_est <- function(x, y, z, s_max, theta_init = NULL, est_method = "gmm") {

    if (est_method == "gmm") {
        ccrm_est_result <- ccrm_est_K2(x, y, z,
                                       theta_init = theta_init, s_max = s_max)
    } else if (est_method == "gmm_sym") {
        ccrm_est_result <- ccrm_est_K2_sym(x, y, z,
                                           theta_init = theta_init, s_max = s_max)
    } else if (est_method == "md") {
        ccrm_est_result <- ccrm_est_hetero(x, y, z,
                                           theta_init = theta_init, s_max = s_max)
    } else {
        stop("Wrong method input.")
    }

    ccrm_est_result
}
