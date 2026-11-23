# Sun-Cai benchmark: unchanged SC functions from the verified local implementation.
sc_eps_estimate <- function(z, mu = 0, sigma = 1) {
  standardized <- (z - mu) / sigma
  xi <- seq(0, 1, by = 0.01)
  t_grid <- seq(0, sqrt(log(length(z))), by = 0.1)

  eps_hat <- vapply(t_grid, function(t) {
    weights <- 1 - abs(xi)
    inflation <- exp((t * xi)^2 / 2)
    char_real <- vapply(xi, function(x) mean(cos(t * x * standardized)), numeric(1))
    1 - sum(weights * inflation * char_real) / sum(weights)
  }, numeric(1))

  eps <- max(eps_hat, na.rm = TRUE)
  min(max(eps, 0), 1)
}

sc_prepare <- function(z, gamma = 0.1, model = c("Gauss01")) {
  model <- match.arg(model)
  if (!length(z) || anyNA(z)) {
    stop("z must be a non-empty vector without missing values.")
  }

  dens_obj <- density(z, from = min(z) - 10, to = max(z) + 10, n = 1000)
  mixture_density <- approx(
    x = dens_obj$x,
    y = dens_obj$y,
    xout = z,
    rule = 2
  )$y

  if (model == "Gauss01") {
    mu <- 0
    sigma <- 1
  }

  pi0 <- 1 - sc_eps_estimate(z, mu = mu, sigma = sigma)
  pi0 <- min(max(pi0, 0), 1)
  lfdr <- pi0 * dnorm(z, mean = mu, sd = sigma) / pmax(mixture_density, .Machine$double.eps)
  lfdr <- pmin(pmax(lfdr, 0), 1)

  list(
    lfdr = lfdr,
    pi0 = pi0,
    mu = mu,
    sigma = sigma,
    gamma = gamma,
    model = model
  )
}

sc_apply <- function(prep, alpha) {
  lfdr <- prep$lfdr
  ord <- order(lfdr)
  lfdr_sorted <- lfdr[ord]
  cumulative_mean <- cumsum(lfdr_sorted) / seq_along(lfdr_sorted)
  k <- max(which(cumulative_mean <= alpha), 0L)

  calls <- rep(FALSE, length(lfdr))
  if (k > 0L) {
    calls[ord[seq_len(k)]] <- TRUE
  }
  calls
}

