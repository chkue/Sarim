#' Generate the structure/penalty matrix for sarim formulae
#' 
#' @description Function used in definition of smooth terms within sarim model formulae.
#' 
#' @param formula usual formula, y ~ x1 + sx(x2), be aware an intercept can only be
#'                included by setting intercept = "TRUE".
#' @param data a data.frame.
#' @param intercept decide if an intercept should be included or not.
#' @param nIter number of MCMC-Iteration.
#' @param burnin number of Burn-In iterations (currently NOT available).
#' @param family family can currently only be "gaussian", "binomial" or "poisson".
#' @param link link-function, currently only be "identity", "logit" or "log" implemented,
#'             if family = "binomial", please choose "logit", if family = "poisson", please choose "log".
#' @param sigma start value for the dispersion parameter, it is assumed that sigma ~ IG(sigma_a, sigma_b),
#'              i.e. it is assumed that sigma is Inverse Gamma distributed with parameters sigma_a and sigma_b.
#' @param sigma_a prior-value, see sigma.
#' @param sigma_b prior-value, see sigma.
#' @param m number of maximal Lanczos-iteration.
#' @param thr threshold, when the Lanczos-algorithm should break.
#' 
#' @export
sarim <- function(formula, data = list(), intercept = "TRUE", nIter = 1000L, burnin = 100L,
                  family = "gaussian", link = "identity",
                  sigma = 1L, sigma_a = 10, sigma_b = 0.1,
                  m = 50L, thr = 0.000001) {
    require(Matrix)
    mf <- stats::model.frame(formula = formula, data = data)
    y <- as.numeric(stats::model.response(mf))
    
    Z <- list()
    K <- list()
    kappaList <- list()
    gammaList <- list()
    solverList <- list()
    kappa_startList <- list()
    etaList <- list()
    
    check_intercept <- list("TRUE", "FALSE")
    if (!is.element(intercept, check_intercept))
        stop('Intercept can only be "TRUE" or "FALSE" ')
    
    check_family <- list("gaussian", "binomial", "poisson")
    if (!is.element(family, check_family))
        stop('"family" can only be "gaussian", "poisson" or "binomial"')
    
    check_link <- list("identity", "log", "logit")
    if (!is.element(link, check_link))
        stop('"link" can only be "identity" (gaussian), "log" (poisson) or "logit" (binomial) ')
    
    
    if ( (family != "gaussian") && (link == "identity") )
        stop('Gibbs sampling only possible if family = "gaussian" ')
    
    
    
    for (i in 2:length(mf)) {
        if (intercept == "TRUE") {
            kappa_startList[[1]] <- as.numeric(1)
            solverList[[1]] <- "rue"
            if (family == "gaussian") {
                gammaList[[1]] <- c(stats::rnorm(1, sd = 0.1))
            }
            if (family != "gaussian") {
                gammaList[[1]] <- c(stats::runif(1, 0, 0.1))
            }
            kappaList[[1]] <- c(1, 0.0001)
            K[[1]] <- as(as.matrix(1), "dgCMatrix")
            Z[[1]] <- as(as.matrix(rep(1, length(y)), nrow = length(y)), "dgCMatrix")
            
            if (!is.vector(mf[[i]])) {
                kappa_startList[[i]] <- attr(mf[[i]], "ka_start")
                solverList[[i]] <- attr(mf[[i]], "solver")
                kappaList[[i]] <- c(attr(mf[[i]], "ka_a"), attr(mf[[i]], "ka_b"))
                
                
                if (is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i]] <- c(stats::rnorm(ncol(mf[[i]]), sd = 0.1))
                    } else {
                        gammaList[[i]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i]] <- as(K_matrix(dimension = dim(mf[[i]])[2],
                                          penalty = attr(mf[[i]], "penalty"),
                                          data = attr(mf[[i]], "K")), "dgCMatrix")
                    attributes(mf[[i]]) <- attributes(mf[[i]])["dim"]
                    Z[[i]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                }
                
                if (!is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i]] <- c( stats::rnorm(ncol(attr(mf[[i]], "Z")), sd = 0.1) )
                    } else {
                        gammaList[[i]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i]] <- as(attr(mf[[i]], "K"), "dgCMatrix")
                    Z[[i]] <- as(attr(mf[[i]], "Z"), "dgCMatrix")
                }
                etaList[[i]] <- Z[[i]] %*% gammaList[[i]]
            }
            if (is.vector(mf[[i]])) {
                kappa_startList[[i]] <- attr(mf[[i]], "ka_start")
                solverList[[i]] <- "rue"
                gammaList[[i]] <- c(stats::rnorm(1, sd = 0.1))
                kappaList[[i]] <- c(1, 0.0001)
                K[[i]] <- as(as.matrix(1), "dgCMatrix")
                Z[[i]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                etaList[[i]] <- Z[[i]] %*% gammaList[[i]]
            }
        } else {
            if (!is.vector(mf[[i]])) {
                kappa_startList[[i - 1]] <- attr(mf[[i]], "ka_start")
                solverList[[i - 1]] <- attr(mf[[i]], "solver")
                gammaList[[i - 1]] <- c(stats::rnorm(ncol(mf[[i]])))
                kappaList[[i - 1]] <- c(attr(mf[[i]], "ka_a"), attr(mf[[i]], "ka_b"))
                
                
                if (is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i - 1]] <- c( stats::rnorm( ncol(mf[[i]]), sd = 0.1 ) )
                    } else {
                        gammaList[[i - 1]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i - 1]] <- as(K_matrix(dimension = dim(mf[[i]])[2],
                                              penalty = attr(mf[[i]], "penalty"),
                                              data = attr(mf[[i]], "K")), "dgCMatrix")
                    
                    attributes(mf[[i]]) <- attributes(mf[[i]])["dim"]
                    Z[[i - 1]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                }
                
                if (!is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i - 1]] <- c( stats::rnorm( ncol(attr( mf[[i]], "Z" )), sd = 0.1 ) )
                    } else {
                        gammaList[[i - 1]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i - 1]] <- as(attr(mf[[i]], "K"), "dgCMatrix")
                    Z[[i - 1]] <- as(attr(mf[[i]], "Z"), "dgCMatrix")
                }
                etaList[[i - 1]] <- Z[[i - 1]] %*% gammaList[[i - 1]]
                
            }
            if (is.vector(mf[[i]])) {
                kappa_startList[[i - 1]] <- as.numeric(1)
                solverList[[i - 1]] <- "rue"
                gammaList[[i - 1]] <- c(stats::rnorm(1, sd = 0.1))
                kappaList[[i - 1]] <- c(1, 0.0001)
                K[[i - 1]] <- as(as.matrix(1), "dgCMatrix")
                Z[[i - 1]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                etaList[[i - 1]] <- Z[[i - 1]] %*% gammaList[[i - 1]]
            }
        }
        
        
    }
    
    sigmavalues <- c(sigma_a, sigma_b)
    sigma <- as.numeric(sigma)
    
    # calculate eta the first time
    eta_first <- rep(0, length(y))
    for (i in 1:length(Z)) eta_first = eta_first + etaList[[i]]
    eta_first <- as.numeric(eta_first)
    
    
    # Using Gibbs or Gibbs-with-MH depending on family
    if (family == "gaussian") {
        out <- Sarim::sarim_gibbs(y, eta_first, Z = Z, K = K, gamma = gammaList,
                                  ka_start = kappa_startList, ka_values = kappaList,
                                  solver = solverList,
                                  sigma = sigma, sigma_values = sigmavalues,
                                  nIter = nIter,
                                  m = m, thr = thr)
        list_out <- list("coef_results" = out$coef_results,
                         "kappa_results" = out$kappa_results,
                         "sigma_results" = out$sigma_results)
    }
    if (family != "gaussian") {
        out <- Sarim::sarim_mcmc(y, eta_first, Z = Z, K = K, gamma = gammaList,
                                 ka_start = kappa_startList, ka_values = kappaList,
                                 solver = solverList, family = family, link = link,
                                 nIter = nIter,
                                 m = m, thr = thr)
        
        list_out <- list("coef_results" = out$coef_results,
                         "kappa_results" = out$kappa_results,
                         "mu_results" = out$mu_results)
    }
    
    
    return(list_out)
}