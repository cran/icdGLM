#' Complete incomplete data
#'
#' This function fills all incomplete data with a set of possible values equally weighted. This is done in order to apply \emph{\link{icdglm}}.
#'
#' @usage expand_data(data, y, missing.x, value.set, weights = rep.int(1, NROW(data)),
#'                    indicator = rep.int(0, NROW(data)))
#' @param data a vector, matrix, list or data frame containing numerics. This data is checked for incompleteness and needs to contain the independent variables for a subsequent regression with n observations and k regressors. Each gap is filled with all values from \code{value.set}. New observations are added for each possible value.
#' @param y a vector of integers or numerics. This vector has to be complete and is the dependent variable for a subsequent regression.
#' @param missing.x a vector that contains integers and gives the position of the independent variables, for which the data will be checked for incompleteness, i.e. for a matrix the position of the corresponding columns.
#' @param value.set a vector of numerics containing all possible values the missing data can take. This set has to be finite.
#' @param weights a vector of numerics giving the initial weight of each observation. Default is 1 for each observation.
#' @param indicator a vector of integers that indicates which observations belong to each other. If some columns with incomplete data were already completed, this vector has to be passed here. For raw incomplete data, the function connects observations which belong to each other. Default is 0 for this vector indicating no connection.
#' @return \code{expand_data} returns a list with the following elements:
#' \itemize{
#' \item{data}{a data frame of the expanded data with all possible observations (independent variables). The dependent variable is included in the last column.}
#' \item{weights}{the weights for each possible observation.}
#' \item{indicator}{a vector which indicates which observations belong to each other. Such observations have the same integer being the indicator.}
#' }
#' @examples data(TLI.data)
#'           expand_data(data = TLI.data[,1:3],
#'           y = TLI.data[,4],
#'           missing.x = 1:3,
#'           value.set = 0:1)
#' @export
expand_data <- function(data,
                        y,
                        missing.x,
                        value.set,
                        weights = rep.int(1, NROW(data)),
                        indicator = rep.int(0, NROW(data))) {

  data <- as.matrix(data)
  y <- as.numeric(y)
  missing.x <- as.numeric(missing.x)
  value.set <- as.numeric(value.set)
  weights <- as.numeric(weights)
  indicator <- as.numeric(indicator)

  complete.data <- data[rowSums(is.na(data)) == 0,]
  complete.y <- y[rowSums(is.na(data)) == 0]
  complete.weights <- weights[rowSums(is.na(data)) == 0]
  complete.indicator <- indicator[rowSums(is.na(data)) == 0]

  incomplete.data <- data[rowSums(is.na(data)) > 0,]
  incomplete.y <- y[rowSums(is.na(data)) > 0]
  incomplete.weights <- weights[rowSums(is.na(data)) > 0]
  incomplete.indicator <- indicator[rowSums(is.na(data)) > 0]

  if (sum(incomplete.indicator) == 0) {
    incomplete.indicator <- as.numeric(rownames(incomplete.data))
  }

  expand.obs.func <- function(row, miss.x, val.set){
    sapply(val.set,
           function(x){
             row[miss.x] <- x
             row[length(row)] <- row[length(row)] / length(val.set)
             return(row)
           })
    }

  if (NROW(incomplete.data) != 0) {
    for (l in seq_along(missing.x)) {
      h <- missing.x[l]

      ## if incom.data is one observation only, incom.data would not be matrix,
      ## therefore define matrix
      #incom.data <- incomplete.data[is.na(incomplete.data[,h])*1 == 1,]
      incom.data <- matrix(incomplete.data[is.na(incomplete.data[, h]) == 1,],
                           ncol = NCOL(data))
      colnames(incom.data) <- colnames(data)

      if (NROW(incom.data) != 0) {
        incom.y <- incomplete.y[is.na(incomplete.data[, h]) == 1]
        incom.weights <- incomplete.weights[is.na(incomplete.data[, h]) == 1]
        incom.indicator <- incomplete.indicator[is.na(incomplete.data[, h]) == 1]

        com.data <- matrix(incomplete.data[is.na(incomplete.data[, h]) != 1,],
                           ncol = NCOL(data))
        colnames(com.data) <- colnames(data)
        #com.data <- incomplete.data[is.na(incomplete.data[,h])*1 != 1,]
        com.y <- incomplete.y[is.na(incomplete.data[, h]) != 1]
        com.weights <- incomplete.weights[is.na(incomplete.data[, h]) != 1]
        com.indicator <- incomplete.indicator[is.na(incomplete.data[, h]) != 1]

        dat <- cbind(incom.data, incom.y, incom.indicator, incom.weights)


        expand <- t(apply(dat, 1, expand.obs.func, val.set = value.set,
                          miss.x = missing.x[l]))

        expand.dat <- apply(array(expand, dim = c(NROW(dat), NCOL(dat),
                                                  length(value.set))), 2, c)

        expand.x <- expand.dat[, 1:NCOL(incom.data)]
        expand.y <- expand.dat[, NCOL(incom.data) + 1]
        expand.indicator <- expand.dat[, NCOL(incom.data) + 2]
        expand.weights <- expand.dat[, NCOL(incom.data) + 3]

        incomplete.data <- rbind(com.data, expand.x)
        incomplete.y <- c(com.y, expand.y)
        incomplete.weights <- c(com.weights, expand.weights)
        incomplete.indicator <- c(com.indicator, expand.indicator)
      }
    }
  }

  if (NROW(incomplete.data) != 0) {
    data <- rbind(incomplete.data, complete.data)
    y <- c(incomplete.y, complete.y)
    weights <- c(incomplete.weights, complete.weights)
    indicator <- c(incomplete.indicator, complete.indicator)

    row.names(data) <- 1:NROW(data)
  }

  data <- cbind(data, y)
  data <- as.data.frame(data)

  expand_data <- list(
    data = data,
    weights = weights,
    indicator = indicator
  )
  return(expand_data)
}
