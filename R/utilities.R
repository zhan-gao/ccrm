initial_value_gap <- function(Eb_1, Eb_2) {

    # Using the estimated variance
    # If it's a valid estimate, return the s.d.
    # If not, return a random number between 0 and 1

    var_hat <- Eb_2 - Eb_1^2
    if (var_hat > 1e-4) {
        return(sqrt(var_hat))
    } else {
        return(runif(1))
    }
}

check_inverse <- function(m){
    "matrix" %in% class(try(solve(m),silent=TRUE))
}
