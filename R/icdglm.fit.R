#' EM by the Method of Weights for Incomplete Data in GLMs (Algorithm)
#'
#' This function applies the EM algorithm by the method of weights to incomplete data in a general linearized model.
#'
#' @usage icdglm.fit(x, y, weights = rep.int(1, NROW(x)), indicator = rep.int(0, NROW(x)),
#'        family = binomial(link = "logit"), control=list())
#' @param x a vector, matrix, list or data frame containing the independent variables
#' @param y a vector of integers or numerics. This is the dependent variable.
#' @param weights a vector which attaches a weight to each observation. For incomplete data, this is obtained from \code{\link{expand_data}}.
#' @param indicator a vector that indicates which observations belong to each other. This is obtained from \code{\link{expand_data}}.
#' @param family family for glm.fit. See \code{\link{glm}} or \code{\link{family}} for more details. Default is \code{binomial(link = "logit")}.
#' @param control a list of control characteristics. See \code{\link{glm.control}} for further information. Default settings are: \code{epsilon = 1e-10, maxit = 100, trace = FALSE}.
#' @return \code{icdglm.fit} returns a list with the following elements:
#' \itemize{
#' \item{x}{a matrix of numerics containing all independent variables}
#' \item{y}{a vector of numerics containing the dependent variable}
#' \item{new.weights}{the new weights obtained in the final iteration of \emph{icdglm.fit}}
#' \item{indicator}{a vector of integers indicating which observations belong to each other}
#' \item{glm.fit.data}{typical \code{glm.fit} output for the last iteration. See \code{\link{glm.fit}} for further information.}
#' \item{coefficients}{a named vector of coefficients}
#' \item{qr}{QR Decomposition of the information matrix}
#' \item{residuals}{the residuals of the final iteration}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#' \item{rank}{the numeric rank of the fitted linear model}
#' \item{family}{the \link{family} object used.}
#' \item{linear.predictors}{the linear fit on link scale}
#' \item{deviance}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero.}
#' \item{aic}{see \link{glm}}
#' \item{null.deviance}{The deviance for the null model, comparable with deviance. The null model will include the offset, and an intercept if there is one in the model. Note that this will be incorrect if the link function depends on the data other than through the fitted mean: specify a zero offset to force a correct calculation.}
#' \item{iter}{an integer containing the number of iterations in \emph{icdglm.fit} before convergence}
#' \item{weights}{the working weights, that is the weights in the final iteration of the IWLS fit.}
#' \item{prior.weights}{the weights initially supplied, a vector of 1s if none were.}
#' \item{df.residual}{the residual degrees of freedom from the initial data set}
#' \item{df.null}{the residual degrees of freedom from initial data set for the null model}
#' \item{model}{model frame}
#' \item{converged}{TRUE if \emph{icdglm} converged.}
#' \item{call}{the match call}
#' \item{formula}{the formula supplied}
#' \item{terms}{the \emph{\link{terms}} object used}
#' \item{data}{the data argument}
#' \item{control}{the value of the \emph{control} argument used}
#' }
#' @references Ibrahim, Joseph G. (1990). \emph{Incomplete Data in Generalized Linear Models}. Journal of the American Statistical Association, Vol.85, No. 411, pp. 765 - 769.
#' @seealso \code{\link{expand_data}}, \code{\link{icdglm}}
#' \code{\link{glm.fit}}, \code{\link{glm.control}}, \code{\link{summary.glm}}
#' @examples data(TLI.data)
#'           complete.data <- expand_data(data = TLI.data[, 1:3],
#'                                        y = TLI.data[, 4],
#'                                        missing.x = 1:3,
#'                                        value.set = 0:1)
#'           example1 <- icdglm.fit(x = complete.data$data[, 1:3],
#'                                  y = complete.data$data[, 4],
#'                                  weights = complete.data$weights,
#'                                  indicator = complete.data$indicator,
#'                                  family = binomial(link = "logit"),
#'                                  control = list(epsilon = 1e-10,
#'                                                 maxit = 100, trace = TRUE))
#' @export
#' @import stats
#' @import Matrix
icdglm.fit <- function(x,
                       y,
                       weights = rep.int(1, NROW(x)),
                       indicator = rep.int(0, NROW(x)),
                       family = binomial(link = "logit"),
                       control = list()) {
  if (any(is.na(x))) {
    stop("icdglm.fit: NAs in x. Please apply expand_data first!", call. = FALSE)
  }

  control <- do.call("glm.control", control)

  x <- as.matrix(x)
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  indicator <- as.numeric(indicator)

  dat <- apply(x, 1, paste, collapse = ";")

  incomplete.obs <- weights != 1

  incomplete.indexes <- which(incomplete.obs)

  sum.weights <- sum(weights)

  conv <- FALSE
  dev.old <- 0

  group.indicator <- rep.int(NA, length(indicator))

  unique.indicator <- unique(indicator)

  for (i in incomplete.indexes) {
    group.indicator[i] <- which(unique.indicator == indicator[i])
  }

  weights.x.obs <- as.numeric(rep.int(NA, NROW(x)))

  matchdat <- match(dat, unique(dat))

  for (iter in 1L:control$maxit) {

    iterations <- iter

    fit1 <- suppressWarnings(glm.fit(x = x, y = y,
                                     family = family, weights = weights))

    if (fit1$converged == FALSE) {
      warning("glm.fit: algorithm did not converge", call. = FALSE)
      break
    }

  if (family$family == "binomial") {
    weights.y <- predict.glm(fit1, type = "response")
    weights.y[y == 0] <- 1 - weights.y[y == 0]
  } else if (family$family == "Gamma") {
    weights.y <- predict.glm(fit1, type = "response")
    print(weights.y)
    weights.y <- dgamma(family$linkinv(y - weights.y), 1) # FIXME
    print(weights.y)
    stop()
  } else {
    # FIXME
  }

    ## very slow for large number of covariates???
    sum1 <- unlist(lapply(unique(dat)[matchdat],
                          function(x)
                          sum(weights[match(dat, x, 0) > 0])))

    weights.x.obs <- sum1 / sum.weights

    product.weights <- weights.x.obs * weights.y

    sum.group.weights <- rowsum(product.weights, group = indicator,
                          reorder = FALSE)

    for (i in incomplete.indexes) {
        sum2 <- sum.group.weights[group.indicator[i],]
        weights[i] <- product.weights[i] / sum2
    }


    dev.new <- fit1$deviance

    dev <- abs(dev.new - dev.old) / (0.1 + abs(dev.new))
    if (control$trace) {
      cat("Deviance deviation = ", dev, " Iterations = ", iter,
          "\n", sep = "")
    }

    if (abs(dev.new - dev.old) / (0.1 + abs(dev.new)) < control$epsilon) {
      conv <- TRUE
      break
    } else {
      dev.old <- dev.new
    }
  }

  if (!conv) {
    warning("icdglm.fit: algorithm did not converge", call. = FALSE)
  }

  Phi <- summary.glm(fit1)$dispersion
  X.tilde <- as.matrix(x)
  W <- Diagonal(length(y), weights)
  W.inv <- Diagonal(length(y), 1 / weights)
  M.tilde <- Diagonal(length(y), fit1$weights)
  M.tilde.2 <- Diagonal(length(y), fit1$weights^2)
  V.tilde.inv.2 <- Diagonal(length(y), (1 / family$variance(fit1$fitted.values))^2)
  H.tilde.2 <- Diagonal(length(y), (y - fit1$fitted.values)^2)

  #info.mat <- t(X.tilde) %*% M.tilde %*% ( Diagonal(length(y),1)  -
  #  ((1/Phi) * M.tilde %*% W.inv %*% V.tilde.inv.2 %*% H.tilde.2 %*%
  #     (Diagonal(length(y),1) - W))) %*% X.tilde

  info.mat <- t(X.tilde) %*% M.tilde %*% (
  Diagonal(length(y), 1 - ((fit1$weights / Phi * (1 - weights) / weights *
                              ((y - fit1$fitted.values) / family$variance(fit1$fitted.values))^2)))) %*% X.tilde

  info.mat <- as.matrix(info.mat)

  print(iter)

  print(family)

  print(y)

  print(weights)

  print(fit1$weights)

  print(fit1$fitted.values)

  print(info.mat)

  qr1 <- list(qr = chol(info.mat),
              pivot = 1:(fit1$rank),
              rank = fit1$rank,
              class = "qr")

  icdglm.fit <- list(x = x,
                     y = fit1$y,
                     new.weights = weights,
                     indicator = indicator,
                     glm.fit.data = fit1,
                     coefficients = fit1$coefficients,
                     qr = qr1,
                     residuals = fit1$residuals,
                     fitted.values = fit1$fitted.values,
                     effects = fit1$effects,
                     R = fit1$R,
                     rank = fit1$rank,
                     family = fit1$family,
                     linear.predictors = fit1$linear.predictors,
                     deviance = fit1$deviance,
                     aic = fit1$aic,
                     null.deviance = fit1$null.deviance,
                     iter = iterations,
                     weights = fit1$weights,
                     prior.weights = fit1$prior.weights,
                     df.residual = sum.weights - NCOL(x),
                     df.null = sum.weights - 1,
                     converged = conv,
                     boundary = fit1$boundary
  )

  return(icdglm.fit)
}
