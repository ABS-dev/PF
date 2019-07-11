
# error checking for use in functions where y can either be a 
# frequency vector of length 4 or a 2x2 matrix.
# checks for type and dim.
.check_y_input_freq <- function(y) {
  if (is.vector(y) & length(y) != 4) {
    stop("'y' input vector incorrect. Must be length 4.")
  } else if (is.matrix(y) & !all(dim(y) == c(2, 2))) {
    stop("'y' matrix must be dim 2 x 2")
  } else if (!(is.vector(y) | is.matrix(y))) {
    stop("'y' must be matrix or vector")
  }
}


# Error checking for use in functions where user can supply frequency vector or
# matrix, or the function calculates it from data+formula.
# Looks for the correct input (either y or data+formula)
# Checks that in the case of y, structure is correct
# This will shortcut any further computation if user did not supply correct info
.check_3input_cases_freq <- function(y = NULL, data = NULL, formula = NULL) {
  
  if (is.null(y) & (is.null(data) | is.null(formula))) {
    # case when y doesn't exist, and either data or formula is missing
    # expect data & formula both present to retrieve values
    stop("Cannot parse input. Are both 'data' and 'formula' present?")
  } else if (!is.null(y) & !(is.null(data) & is.null(formula))) {
    # case when y is present, both data and formula must be null
    stop("Cannot parse input. If 'y' is included, both 'formula' and 'data' ",
      "must be NULL.")
  } else if (!is.null(y)) {
    # case when only the y value is supplied
    .check_y_input_freq(y)
  }
  
}
