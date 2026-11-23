fit_noisysbmGGM_once <- function(X,
                                 NIG = FALSE,
                                 threshold = 0.5,
                                 Nbrepet = noisysbmGGM_Nbrepet,
                                 rho = NULL,
                                 tau = NULL,
                                 a = NULL,
                                 b = NULL,
                                 c = NULL,
                                 d = NULL,
                                 n0 = 1,
                                 eta0 = 1,
                                 zeta0 = 1,
                                 Qup = noisysbmGGM_Qup,
                                 nbCores = noisysbmGGM_nbCores,
                                 nbOfZ = noisysbmGGM_nbOfZ,
                                 sigma0 = 1,
                                 sigma1 = 1,
                                 percentageOfPerturbation = 0.3,
                                 verbatim = FALSE) {
  p <- ncol(X)
  if (nrow(X) != p) {
    stop("X must be a square matrix.")
  }

  dataVec <- X[lower.tri(X)]

  if (NIG) {
    if (is.null(a)) a <- 0
    if (is.null(b)) b <- 1
    if (is.null(c)) c <- 1
    if (is.null(d)) d <- 1

    fit <- mainSearchOpti_NIG(
      X, Qup = Qup, threshold = threshold, Nbrepet = Nbrepet,
      nbCores = nbCores, nbOfZ = nbOfZ,
      percentageOfPerturbation = percentageOfPerturbation,
      a = a, b = b, c = c, d = d, n0 = n0, eta0 = eta0, zeta0 = zeta0,
      sigma0 = sigma0, fast = TRUE, verbatim = verbatim
    )
    selected_icl <- fit$ICL
    Q_before_merge <- fit$Q

    if (fit$Q > 1) {
      merged <- Merge_Nondirige_NIG(
        dataVec, fit$Z, fit$Zmatrix, fit$Q, fit$Rho, fit$A, fit$theta,
        threshold = threshold, a = a, b = b, c = c, d = d,
        n0 = n0, eta0 = eta0, zeta0 = zeta0
      )
      fit$Rho <- merged$Rhomerge
      fit$Q <- merged$Qmerge
      fit$Z <- merged$Zmerge
      fit$theta <- merged$thetamerge
    }
  } else {
    if (is.null(rho)) rho <- 1
    if (is.null(tau)) tau <- 1

    fit <- mainSearchOpti(
      X, Qup = Qup, threshold = threshold, Nbrepet = Nbrepet,
      nbCores = nbCores, nbOfZ = nbOfZ,
      percentageOfPerturbation = percentageOfPerturbation,
      rho = rho, tau = tau, n0 = n0, eta0 = eta0, zeta0 = zeta0,
      sigma0 = sigma0, sigma1 = sigma1, fast = TRUE, verbatim = verbatim
    )
    selected_icl <- fit$ICL
    Q_before_merge <- fit$Q

    if (fit$Q > 1) {
      merged <- Merge_Nondirige(
        dataVec, fit$Z, fit$Zmatrix, fit$Q, fit$Rho, fit$A, fit$theta,
        threshold = threshold, rho = rho, tau = tau,
        n0 = n0, eta0 = eta0, zeta0 = zeta0
      )
      fit$Rho <- merged$Rhomerge
      fit$Q <- merged$Qmerge
      fit$Z <- merged$Zmerge
      fit$theta <- merged$thetamerge
    }
  }

  list(
    Rho = fit$Rho, Z = fit$Z, theta = fit$theta, Q = fit$Q,
    ICL = selected_icl, Q_before_merge = Q_before_merge, dataVec = dataVec
  )
}

infer_noisysbmGGM_graph <- function(fit, alpha) {
  InferGraph(fit$Rho, fit$theta, fit$Z, fit$dataVec, alpha)
}
