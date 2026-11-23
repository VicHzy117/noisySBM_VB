#ifndef NSBM_VI_VARIATIONAL_BLOCKS_H
#define NSBM_VI_VARIATIONAL_BLOCKS_H

#include <Rcpp.h>

namespace nsbmvi {

inline int block_count(const int K) {
  if (K < 1) {
    Rcpp::stop("K must be positive.");
  }
  return K * (K + 1) / 2;
}

inline int block_index(int r, int s, const int K) {
  if (r > s) {
    std::swap(r, s);
  }
  if (r < 0 || s >= K) {
    Rcpp::stop("Block index is out of range.");
  }
  return r * K - r * (r - 1) / 2 + (s - r);
}

inline double stable_logistic(const double value) {
  if (value >= 0.0) {
    return 1.0 / (1.0 + std::exp(-value));
  }
  const double exp_value = std::exp(value);
  return exp_value / (1.0 + exp_value);
}

inline double bernoulli_entropy(const double probability) {
  if (probability <= 0.0 || probability >= 1.0) {
    return 0.0;
  }
  return -probability * std::log(probability) -
    (1.0 - probability) * std::log(1.0 - probability);
}

inline void validate_tau(const Rcpp::NumericMatrix& tau) {
  if (tau.nrow() < 1 || tau.ncol() < 1) {
    Rcpp::stop("tau must have at least one row and one column.");
  }

  for (int i = 0; i < tau.nrow(); ++i) {
    double row_sum = 0.0;
    for (int k = 0; k < tau.ncol(); ++k) {
      const double value = tau(i, k);
      if (!R_finite(value) || value < 0.0 || value > 1.0) {
        Rcpp::stop("Every tau entry must be finite and lie in [0, 1].");
      }
      row_sum += value;
    }
    if (std::abs(row_sum - 1.0) > 1e-8) {
      Rcpp::stop("Every row of tau must sum to one.");
    }
  }
}

inline void validate_X(const Rcpp::NumericMatrix& X) {
  if (X.nrow() < 2 || X.ncol() != X.nrow()) {
    Rcpp::stop("X must be a square matrix with at least two rows.");
  }
  for (int i = 0; i < X.nrow(); ++i) {
    for (int j = i + 1; j < X.ncol(); ++j) {
      if (!R_finite(X(i, j))) {
        Rcpp::stop("The upper triangle of X must be finite.");
      }
    }
  }
}

inline void validate_block_vector(const Rcpp::NumericVector& values,
                                  const int expected_size,
                                  const char* name,
                                  const bool positive) {
  if (values.size() != expected_size) {
    Rcpp::stop("%s has the wrong length.", name);
  }
  for (int b = 0; b < expected_size; ++b) {
    if (!R_finite(values[b]) || (positive && values[b] <= 0.0)) {
      Rcpp::stop("%s contains an invalid value.", name);
    }
  }
}

inline Rcpp::NumericMatrix compute_block_weights(
    const Rcpp::NumericMatrix& tau) {
  validate_tau(tau);
  const int n = tau.nrow();
  const int K = tau.ncol();
  const int B = block_count(K);
  const int M = n * (n - 1) / 2;
  Rcpp::NumericMatrix weights(M, B);

  int edge = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      for (int r = 0; r < K; ++r) {
        for (int s = r; s < K; ++s) {
          const int b = block_index(r, s, K);
          weights(edge, b) = (r == s)
            ? tau(i, r) * tau(j, r)
            : tau(i, r) * tau(j, s) + tau(i, s) * tau(j, r);
        }
      }
      ++edge;
    }
  }

  return weights;
}

inline Rcpp::NumericMatrix update_phi(
    const Rcpp::NumericMatrix& X,
    const int K,
    const Rcpp::NumericVector& beta_a,
    const Rcpp::NumericVector& beta_b,
    const Rcpp::NumericVector& nig_m,
    const Rcpp::NumericVector& nig_kappa,
    const Rcpp::NumericVector& nig_shape,
    const Rcpp::NumericVector& nig_scale) {
  validate_X(X);
  const int B = block_count(K);
  validate_block_vector(beta_a, B, "beta_a", true);
  validate_block_vector(beta_b, B, "beta_b", true);
  validate_block_vector(nig_m, B, "nig_m", false);
  validate_block_vector(nig_kappa, B, "nig_kappa", true);
  validate_block_vector(nig_shape, B, "nig_shape", true);
  validate_block_vector(nig_scale, B, "nig_scale", true);

  Rcpp::NumericVector log_q(B);
  Rcpp::NumericVector log_one_minus_q(B);
  Rcpp::NumericVector expected_log_variance(B);

  for (int b = 0; b < B; ++b) {
    log_q[b] = R::digamma(beta_a[b]) - R::digamma(beta_a[b] + beta_b[b]);
    log_one_minus_q[b] =
      R::digamma(beta_b[b]) - R::digamma(beta_a[b] + beta_b[b]);
    expected_log_variance[b] =
      std::log(nig_scale[b]) - R::digamma(nig_shape[b]);
  }

  const int n = X.nrow();
  const int M = n * (n - 1) / 2;
  Rcpp::NumericMatrix phi(M, B);
  int edge = 0;

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double x = X(i, j);
      const double h0 = R::dnorm(x, 0.0, 1.0, true);

      for (int b = 0; b < B; ++b) {
        const double expected_scaled_square =
          (nig_shape[b] / nig_scale[b]) *
            (x - nig_m[b]) * (x - nig_m[b]) +
          1.0 / nig_kappa[b];
        const double h1 = -0.5 *
          (std::log(2.0 * M_PI) + expected_log_variance[b] +
           expected_scaled_square);
        const double logit =
          log_q[b] - log_one_minus_q[b] + h1 - h0;
        phi(edge, b) = stable_logistic(logit);
      }
      ++edge;
    }
  }

  return phi;
}

inline Rcpp::List update_block_posteriors(
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericMatrix& tau,
    const Rcpp::NumericMatrix& phi,
    const double alpha,
    const double beta,
    const double alpha_prime,
    const double beta_prime,
    const double prior_mean,
    const double prior_kappa,
    const double prior_shape,
    const double prior_scale) {
  validate_X(X);
  validate_tau(tau);
  const int n = X.nrow();
  const int K = tau.ncol();
  const int B = block_count(K);
  const int M = n * (n - 1) / 2;

  if (tau.nrow() != n) {
    Rcpp::stop("tau and X have incompatible dimensions.");
  }
  if (phi.nrow() != M || phi.ncol() != B) {
    Rcpp::stop("phi must have M rows and K(K+1)/2 columns.");
  }
  if (alpha <= 0.0 || beta <= 0.0 ||
      alpha_prime <= 0.0 || beta_prime <= 0.0 ||
      prior_kappa <= 0.0 || prior_shape <= 0.0 || prior_scale <= 0.0) {
    Rcpp::stop("All Beta and NIG scale/shape hyperparameters must be positive.");
  }

  Rcpp::NumericVector active(B);
  Rcpp::NumericVector inactive(B);
  Rcpp::NumericVector sum_x(B);
  Rcpp::NumericVector sum_x2(B);
  int edge = 0;

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double x = X(i, j);
      for (int r = 0; r < K; ++r) {
        for (int s = r; s < K; ++s) {
          const int b = block_index(r, s, K);
          const double weight = (r == s)
            ? tau(i, r) * tau(j, r)
            : tau(i, r) * tau(j, s) + tau(i, s) * tau(j, r);
          const double edge_probability = phi(edge, b);
          if (!R_finite(edge_probability) ||
              edge_probability < 0.0 || edge_probability > 1.0) {
            Rcpp::stop("phi contains an invalid probability.");
          }
          const double active_weight = weight * edge_probability;
          active[b] += active_weight;
          inactive[b] += weight * (1.0 - edge_probability);
          sum_x[b] += active_weight * x;
          sum_x2[b] += active_weight * x * x;
        }
      }
      ++edge;
    }
  }

  Rcpp::NumericVector beta_a(B);
  Rcpp::NumericVector beta_b(B);
  Rcpp::NumericVector nig_m(B);
  Rcpp::NumericVector nig_kappa(B);
  Rcpp::NumericVector nig_shape(B);
  Rcpp::NumericVector nig_scale(B);

  for (int r = 0; r < K; ++r) {
    for (int s = r; s < K; ++s) {
      const int b = block_index(r, s, K);
      const double a0 = (r == s) ? alpha_prime : alpha;
      const double b0 = (r == s) ? beta_prime : beta;
      beta_a[b] = a0 + active[b];
      beta_b[b] = b0 + inactive[b];

      nig_kappa[b] = prior_kappa + active[b];
      nig_m[b] =
        (prior_kappa * prior_mean + sum_x[b]) / nig_kappa[b];
      nig_shape[b] = prior_shape + active[b] / 2.0;
      nig_scale[b] = prior_scale + 0.5 *
        (sum_x2[b] + prior_kappa * prior_mean * prior_mean -
         nig_kappa[b] * nig_m[b] * nig_m[b]);

      if (nig_scale[b] <= 0.0 && nig_scale[b] > -1e-10) {
        nig_scale[b] = std::numeric_limits<double>::min();
      }
      if (!R_finite(nig_scale[b]) || nig_scale[b] <= 0.0) {
        Rcpp::stop("Updated NIG scale is non-positive or non-finite.");
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("beta_a") = beta_a,
    Rcpp::Named("beta_b") = beta_b,
    Rcpp::Named("nig_m") = nig_m,
    Rcpp::Named("nig_kappa") = nig_kappa,
    Rcpp::Named("nig_shape") = nig_shape,
    Rcpp::Named("nig_scale") = nig_scale,
    Rcpp::Named("active") = active,
    Rcpp::Named("inactive") = inactive,
    Rcpp::Named("sum_x") = sum_x,
    Rcpp::Named("sum_x2") = sum_x2
  );
}

inline Rcpp::NumericMatrix compute_edge_scores(
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericMatrix& phi,
    const int K,
    const Rcpp::NumericVector& beta_a,
    const Rcpp::NumericVector& beta_b,
    const Rcpp::NumericVector& nig_m,
    const Rcpp::NumericVector& nig_kappa,
    const Rcpp::NumericVector& nig_shape,
    const Rcpp::NumericVector& nig_scale) {
  validate_X(X);
  const int B = block_count(K);
  const int n = X.nrow();
  const int M = n * (n - 1) / 2;
  if (phi.nrow() != M || phi.ncol() != B) {
    Rcpp::stop("phi has incompatible dimensions.");
  }
  validate_block_vector(beta_a, B, "beta_a", true);
  validate_block_vector(beta_b, B, "beta_b", true);
  validate_block_vector(nig_m, B, "nig_m", false);
  validate_block_vector(nig_kappa, B, "nig_kappa", true);
  validate_block_vector(nig_shape, B, "nig_shape", true);
  validate_block_vector(nig_scale, B, "nig_scale", true);

  Rcpp::NumericMatrix scores(M, B);
  int edge = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double x = X(i, j);
      const double h0 = R::dnorm(x, 0.0, 1.0, true);
      for (int b = 0; b < B; ++b) {
        const double log_q =
          R::digamma(beta_a[b]) - R::digamma(beta_a[b] + beta_b[b]);
        const double log_one_minus_q =
          R::digamma(beta_b[b]) - R::digamma(beta_a[b] + beta_b[b]);
        const double h1 = -0.5 * (
          std::log(2.0 * M_PI) +
          std::log(nig_scale[b]) - R::digamma(nig_shape[b]) +
          (nig_shape[b] / nig_scale[b]) *
            (x - nig_m[b]) * (x - nig_m[b]) +
          1.0 / nig_kappa[b]
        );
        const double p = phi(edge, b);
        scores(edge, b) =
          p * (log_q + h1) +
          (1.0 - p) * (log_one_minus_q + h0) +
          bernoulli_entropy(p);
      }
      ++edge;
    }
  }

  return scores;
}

}  // namespace nsbmvi

#endif

