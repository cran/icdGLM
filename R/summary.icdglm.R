#' Summarizing Output of an EM Algorithm by the Method of Weights Using GLMs
#'
#' This function gives a summary of the output of \emph{\link{icdglm}}. \emph{summary.icdglm} inherits from \emph{summary.glm}.
#'
#' @usage \method{summary}{icdglm}(object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...)
#' @param object an object of class "icdglm", usually, a result of a call to \emph{\link{icdglm}}.
#' @param dispersion the dispersion parameter for the family used. Either a single numerical value or NULL (the default), when it is inferred from object (see details of \emph{\link{summary.glm}}).
#' @param correlation logical, if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor logical, if TRUE, print the correlations in a symbolic form (see \emph{\link{symnum}}) rather than as numbers.
#' @param \dots further arguments passed to or from other methods.
#' @note The description of this function is taken from \emph{\link{summary.glm}} apart from a few differences.
#' @return \emph{summary.icdglm} returns an object of class "\emph{summary.icdglm}", a list with components:
#' \itemize{
#' \item{call}{function call of \emph{object}}
#' \item{terms}{the \emph{\link{terms}} object used.}
#' \item{family}{the component from \emph{object}}
#' \item{deviance}{the component from \emph{object}}
#' \item{aic}{the component from \emph{object}}
#' \item{df.residual}{the residual degrees of freedom of the initial data set}
#' \item{null.deviance}{the component from \emph{object}}
#' \item{df.null}{the residual degrees of freedom for the null model.}
#' \item{iter}{the number of iterations in \emph{icdglm.fit}, component from \emph{object}}
#' \item{deviance.resid}{the deviance residuals: see \emph{\link{residuals.glm}}}
#' \item{coefficients}{the matrix of coefficients, (corrected) standard errors, t-values and p-values.}
#' \item{aliased}{named logical vector showing if the original coefficients are aliased.}
#' \item{dispersion}{either the supplied argument or the inferred/estimated dispersion if the latter is NULL.}
#' \item{df}{a 3-vector of the rank of the model and the number of residual degrees of freedom, plus number of coefficients (including aliased ones).}
#' \item{cov.unscaled}{the unscaled (dispersion = 1) estimated covariance matrix of the estimated coefficients.}
#' \item{cov.scaled}{ditto, scaled by dispersion}
#' \item{correlation}{(only if \emph{correlation} is \emph{TRUE}) The estimated correlations of the estimated coefficients.}
#' \item{symbolic.cor}{(only if \emph{correlation} is \emph{TRUE}) The value of the argument symbolic.cor.}
#' }
#' @seealso \code{\link{icdglm}}, \code{\link{summary.glm}}, \code{\link{summary}}, \code{\link{glm}}
#' @export
#' @import stats
summary.icdglm <-
  function(object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
            ...) {

  ans <- summary.glm(object, dispersion = dispersion,
                     correlation = correlation, symbolic.cor = symbolic.cor,
                     ...)

    class(ans) <- c("summary.icdglm", "summary.glm")
    return(ans)
  }
