# Thin R wrapper around the Rcpp pooled-posterior spectral initializer.
#
# The complete initialization path -- pooled EM, soft adjacency construction,
# eigendecomposition, row normalization, multistart k-means, label softening,
# perturbations, and phi initialization -- lives in
# code/spectral_initialization.cpp.  R retains only the short ELBO screening
# and continuation logic for the fixed-K VI fit.

make_nsbm_spectral_initialization <- function(
    X, K, seed = 1L, epsilon = 0.1, kmeans_nstart = 25L,
    n_perturb = 2L, perturb_fraction = 0.15,
    pooled_control = list()) {
  if (!exists("nsbm_spectral_init_cpp", mode = "function")) {
    stop("Compile spectral_initialization.cpp before calling this initializer")
  }

  controls <- utils::modifyList(
    list(
      null_sd = 1,
      n_starts = 5L,
      max_iter = 5000L,
      tol = 1e-8,
      min_eff_n = 50,
      sigma2_floor = 1e-3,
      probability_floor = 1e-8
    ),
    pooled_control
  )

  nsbm_spectral_init_cpp(
    X = X,
    K = as.integer(K),
    seed = as.integer(seed),
    epsilon = epsilon,
    kmeans_nstart = as.integer(kmeans_nstart),
    n_perturb = as.integer(n_perturb),
    perturb_fraction = perturb_fraction,
    null_sd = controls$null_sd,
    em_n_starts = as.integer(controls$n_starts),
    em_max_iter = as.integer(controls$max_iter),
    em_tolerance = controls$tol,
    min_eff_n = controls$min_eff_n,
    sigma2_floor = controls$sigma2_floor,
    probability_floor = controls$probability_floor
  )
}

fit_nsbm_vi_spectral <- function(
    X, K, hyperparameters, seed = 1L,
    epsilon = 0.1, kmeans_nstart = 25L,
    n_perturb = 2L, perturb_fraction = 0.15,
    screen_iter = 15L, max_iter = 300L,
    tolerance = 1e-6, min_iter = 3L,
    pooled_control = list(), verbose = FALSE, count_pmf_rebuild = TRUE) {
  if (!exists("nsbm_vi_fixed_k_cpp", mode = "function")) {
    stop("Compile nsbm_vi.cpp before calling this wrapper")
  }

  initialization <- make_nsbm_spectral_initialization(
    X = X, K = K, seed = seed, epsilon = epsilon,
    kmeans_nstart = kmeans_nstart,
    n_perturb = n_perturb, perturb_fraction = perturb_fraction,
    pooled_control = pooled_control
  )

  screen_started <- proc.time()
  screen_fits <- lapply(initialization$starts, function(start) {
    nsbm_vi_fixed_k_cpp(
      X, K = K, tau_init = start$tau,
      hyperparameters = hyperparameters,
      max_iter = screen_iter, tolerance = tolerance,
      min_iter = min(min_iter, screen_iter), verbose = FALSE,
      phi_init = initialization$phi, count_pmf_rebuild = count_pmf_rebuild
    )
  })
  screen_elapsed <- unname((proc.time() - screen_started)["elapsed"])
  final_elbos <- vapply(
    screen_fits, function(fit) tail(fit$elbo_trace, 1L), numeric(1L)
  )
  best <- which.max(final_elbos)

  # Continue from the winning screened factors.  The C++ driver reconstructs
  # the conjugate block factors from this tau/phi pair before resuming CAVI.
  fit_started <- proc.time()
  fit <- nsbm_vi_fixed_k_cpp(
    X, K = K,
    tau_init = screen_fits[[best]]$tau,
    hyperparameters = hyperparameters,
    max_iter = max_iter, tolerance = tolerance,
    min_iter = min_iter, verbose = verbose,
    phi_init = screen_fits[[best]]$phi, count_pmf_rebuild = count_pmf_rebuild
  )
  fit_elapsed <- unname((proc.time() - fit_started)["elapsed"])

  screening <- data.frame(
    start = names(initialization$starts),
    iterations = vapply(screen_fits, `[[`, integer(1L), "iterations"),
    final_elbo = final_elbos,
    selected = seq_along(screen_fits) == best,
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    initialization = initialization,
    screening = screening,
    selected_start = names(initialization$starts)[best],
    selected_initial_labels = initialization$starts[[best]]$labels,
    timing = c(
      initialization = initialization$timing[["total"]],
      screening = screen_elapsed,
      final_fit = fit_elapsed,
      total = initialization$timing[["total"]] + screen_elapsed + fit_elapsed
    )
  )
}

# Fit a path of fixed-K variational models and select K by the maximized ELBO:
# score(K) = ELBO(K).  log_p_K is retained in the returned path only as a
# diagnostic so old ELBO + log p_K(K) selections can still be reconstructed.
fit_nsbm_vi_select_k <- function(
    X, K_candidates = 1:6, hyperparameters,
    lambda_mfm = NULL,
    seed = 1L, epsilon = 0.1, kmeans_nstart = 25L,
    n_perturb = 2L, perturb_fraction = 0.15,
    screen_iter = 15L, max_iter = 300L,
    tolerance = 1e-6, min_iter = 3L,
    pooled_control = list(), verbose = FALSE, count_pmf_rebuild = TRUE) {
  K_candidates <- sort(unique(as.integer(K_candidates)))
  if (!length(K_candidates) || anyNA(K_candidates) ||
      any(K_candidates < 1L) || any(K_candidates > nrow(X))) {
    stop("K_candidates must be distinct integers between 1 and nrow(X)")
  }
  if (!is.null(lambda_mfm) &&
      (!is.finite(lambda_mfm) || lambda_mfm <= 0)) {
    stop("lambda_mfm must be NULL or positive")
  }
  log_p_K_diagnostic <- function(K) {
    if (is.null(lambda_mfm)) NA_real_ else {
      stats::dpois(K - 1L, lambda_mfm, log = TRUE)
    }
  }

  fits <- vector("list", length(K_candidates))
  path_rows <- vector("list", length(K_candidates))
  names(fits) <- as.character(K_candidates)

  path_started <- proc.time()
  for (index in seq_along(K_candidates)) {
    K <- K_candidates[[index]]
    candidate_started <- proc.time()
    candidate <- tryCatch(
      fit_nsbm_vi_spectral(
        X = X, K = K, hyperparameters = hyperparameters,
        seed = as.integer(seed + 1009L * K),
        epsilon = epsilon, kmeans_nstart = kmeans_nstart,
        n_perturb = if (K == 1L) 0L else n_perturb,
        perturb_fraction = perturb_fraction,
        screen_iter = screen_iter, max_iter = max_iter,
        tolerance = tolerance, min_iter = min_iter,
        pooled_control = pooled_control, verbose = verbose, count_pmf_rebuild = count_pmf_rebuild
      ),
      error = function(error) error
    )
    elapsed <- unname((proc.time() - candidate_started)["elapsed"])

    if (inherits(candidate, "error")) {
      path_rows[[index]] <- data.frame(
        K = K, success = FALSE, error = conditionMessage(candidate),
        converged = FALSE, monotone = FALSE, iterations = NA_integer_,
        selected_start = NA_character_, final_elbo = NA_real_,
        log_p_K = log_p_K_diagnostic(K),
        score = -Inf, model_weight = 0,
        initialization_sec = NA_real_, screening_sec = NA_real_,
        final_fit_sec = NA_real_, total_method_sec = elapsed,
        stringsAsFactors = FALSE
      )
      next
    }

    fits[[index]] <- candidate
    final_elbo <- tail(candidate$fit$elbo_trace, 1L)
    log_p_K <- log_p_K_diagnostic(K)
    path_rows[[index]] <- data.frame(
      K = K, success = TRUE, error = NA_character_,
      converged = candidate$fit$converged,
      monotone = candidate$fit$monotone,
      iterations = candidate$fit$iterations,
      selected_start = candidate$selected_start,
      final_elbo = final_elbo,
      log_p_K = log_p_K,
      score = final_elbo,
      model_weight = NA_real_,
      initialization_sec = candidate$timing[["initialization"]],
      screening_sec = candidate$timing[["screening"]],
      final_fit_sec = candidate$timing[["final_fit"]],
      total_method_sec = candidate$timing[["total"]],
      stringsAsFactors = FALSE
    )
  }

  path <- do.call(rbind, path_rows)
  valid <- which(path$success & is.finite(path$score))
  if (!length(valid)) stop("Every candidate-K variational fit failed")
  shifted <- path$score[valid] - max(path$score[valid])
  path$model_weight[valid] <- exp(shifted) / sum(exp(shifted))
  path$model_weight[!path$success] <- 0
  selected_index <- valid[[which.max(path$score[valid])]]
  path$selected <- seq_len(nrow(path)) == selected_index

  list(
    selected_K = path$K[[selected_index]],
    selected = fits[[selected_index]],
    path = path,
    fits = fits,
    settings = list(
      K_candidates = K_candidates,
      selection_criterion = "ELBO"
    ),
    timing = c(
      candidate_fits = sum(path$total_method_sec, na.rm = TRUE),
      path_wall = unname((proc.time() - path_started)["elapsed"])
    )
  )
}
