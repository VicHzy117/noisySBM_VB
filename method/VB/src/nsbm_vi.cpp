#include "exact_count.h"
#include "variational_blocks.h"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

struct ElboResult {
  double allocation;
  double edge;
  double beta;
  double nig;

  double total() const {
    return allocation + edge + beta + nig;
  }
};

using LongPmf = std::vector<long double>;

void normalize_long_pmf(LongPmf& pmf) {
  long double total = 0.0L;
  for (std::size_t index = 0; index < pmf.size(); ++index) {
    if (!std::isfinite(pmf[index])) {
      stop("Long-double count pmf contains a non-finite value.");
    }
    if (pmf[index] < -1e-12L) {
      stop(
        "Long-double count pmf contains a materially negative value: %.17g.",
        static_cast<double>(pmf[index])
      );
    }
    if (pmf[index] < 0.0L) {
      pmf[index] = 0.0L;
    }
    total += pmf[index];
  }
  if (!std::isfinite(total) || total <= 0.0L) {
    stop("Long-double count pmf has invalid total mass.");
  }
  for (std::size_t index = 0; index < pmf.size(); ++index) {
    pmf[index] /= total;
  }
}

LongPmf count_pmf_long(const NumericMatrix& tau, const int component,
                       const int excluded_node = -1) {
  const int n = tau.nrow();
  const int size = n + 1 - (excluded_node >= 0 ? 1 : 0);
  LongPmf pmf(size, 0.0L);
  pmf[0] = 1.0L;
  int current_degree = 0;

  for (int i = 0; i < n; ++i) {
    if (i == excluded_node) continue;
    const long double p = static_cast<long double>(tau(i, component));
    LongPmf next(size, 0.0L);
    for (int count = 0; count <= current_degree; ++count) {
      next[count] += (1.0L - p) * pmf[count];
      next[count + 1] += p * pmf[count];
    }
    pmf.swap(next);
    ++current_degree;
  }
  normalize_long_pmf(pmf);
  return pmf;
}

LongPmf remove_factor_long(const LongPmf& full,
                           const double probability) {
  const int full_size = full.size();
  const int reduced_size = full_size - 1;
  const long double p = static_cast<long double>(probability);
  LongPmf reduced(reduced_size, 0.0L);

  if (p <= 0.5L) {
    const long double denominator = 1.0L - p;
    reduced[0] = full[0] / denominator;
    for (int count = 1; count < reduced_size; ++count) {
      reduced[count] =
        (full[count] - p * reduced[count - 1]) / denominator;
    }
  } else {
    reduced[reduced_size - 1] = full[full_size - 1] / p;
    for (int count = reduced_size - 1; count >= 1; --count) {
      reduced[count - 1] =
        (full[count] - (1.0L - p) * reduced[count]) / p;
    }
  }

  normalize_long_pmf(reduced);
  return reduced;
}

LongPmf remove_factor_or_rebuild(const LongPmf& full, const NumericMatrix& tau,
                                 const int component, const int node,
                                 int& numerical_rebuilds) {
  try {
    return remove_factor_long(full, tau(node, component));
  } catch (const std::exception&) {
    // Polynomial division can amplify round-off even with long doubles.
    // Recompute the same leave-one-out distribution by positive-only
    // convolution from current tau; retain the existing validity tolerance.
    // If the source probabilities themselves are invalid, this still fails.
    LongPmf rebuilt = count_pmf_long(tau, component, node);
    ++numerical_rebuilds;
    return rebuilt;
  }
}

LongPmf add_factor_long(const LongPmf& reduced,
                        const double probability) {
  const long double p = static_cast<long double>(probability);
  LongPmf full(reduced.size() + 1, 0.0L);
  for (std::size_t count = 0; count < reduced.size(); ++count) {
    full[count] += (1.0L - p) * reduced[count];
    full[count + 1] += p * reduced[count];
  }
  normalize_long_pmf(full);
  return full;
}

double prior_score_long(const LongPmf& reduced, const double gamma) {
  long double score = 0.0L;
  for (std::size_t count = 0; count < reduced.size(); ++count) {
    score += reduced[count] *
      std::log(static_cast<long double>(count) + gamma);
  }
  return static_cast<double>(score);
}

double get_hyperparameter(const List& hyperparameters,
                          const char* name,
                          const double default_value) {
  if (hyperparameters.containsElementNamed(name)) {
    const double value = as<double>(hyperparameters[name]);
    if (!R_finite(value)) {
      stop("Hyperparameter %s must be finite.", name);
    }
    return value;
  }
  return default_value;
}

NumericMatrix initialize_tau(const int n, const int K) {
  NumericMatrix tau(n, K);
  if (K == 1) {
    std::fill(tau.begin(), tau.end(), 1.0);
    return tau;
  }

  const double main_probability = 0.90;
  const double other_probability =
    (1.0 - main_probability) / static_cast<double>(K - 1);
  for (int i = 0; i < n; ++i) {
    const int label = i % K;
    for (int k = 0; k < K; ++k) {
      tau(i, k) = (k == label) ? main_probability : other_probability;
    }
  }
  return tau;
}

NumericMatrix count_pmf_matrix(const NumericMatrix& tau) {
  const int n = tau.nrow();
  const int K = tau.ncol();
  NumericMatrix output(K, n + 1);

  for (int k = 0; k < K; ++k) {
    NumericVector probabilities(n);
    for (int i = 0; i < n; ++i) {
      probabilities[i] = tau(i, k);
    }
    const NumericVector pmf = nsbmvi::count_pmf(probabilities);
    for (int count = 0; count <= n; ++count) {
      output(k, count) = pmf[count];
    }
  }
  return output;
}

double allocation_elbo(const NumericMatrix& tau,
                       const double gamma) {
  nsbmvi::validate_tau(tau);
  const int n = tau.nrow();
  const int K = tau.ncol();
  const NumericMatrix pmfs = count_pmf_matrix(tau);

  double value =
    R::lgammafn(K * gamma) -
    R::lgammafn(n + K * gamma) -
    K * R::lgammafn(gamma);

  for (int k = 0; k < K; ++k) {
    for (int count = 0; count <= n; ++count) {
      value += pmfs(k, count) *
        R::lgammafn(static_cast<double>(count) + gamma);
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      if (tau(i, k) > 0.0) {
        value -= tau(i, k) * std::log(tau(i, k));
      }
    }
  }

  return value;
}

double edge_elbo(const NumericMatrix& tau,
                 const NumericMatrix& edge_scores) {
  const NumericMatrix weights = nsbmvi::compute_block_weights(tau);
  if (weights.nrow() != edge_scores.nrow() ||
      weights.ncol() != edge_scores.ncol()) {
    stop("edge_scores has incompatible dimensions.");
  }

  double value = 0.0;
  for (int edge = 0; edge < weights.nrow(); ++edge) {
    for (int b = 0; b < weights.ncol(); ++b) {
      value += weights(edge, b) * edge_scores(edge, b);
    }
  }
  return value;
}

double beta_elbo(const int K,
                 const NumericVector& beta_a,
                 const NumericVector& beta_b,
                 const double alpha,
                 const double beta,
                 const double alpha_prime,
                 const double beta_prime) {
  double value = 0.0;
  for (int r = 0; r < K; ++r) {
    for (int s = r; s < K; ++s) {
      const int b = nsbmvi::block_index(r, s, K);
      const double a0 = (r == s) ? alpha_prime : alpha;
      const double b0 = (r == s) ? beta_prime : beta;
      const double u1 =
        R::digamma(beta_a[b]) - R::digamma(beta_a[b] + beta_b[b]);
      const double u0 =
        R::digamma(beta_b[b]) - R::digamma(beta_a[b] + beta_b[b]);
      value +=
        -R::lbeta(a0, b0) + R::lbeta(beta_a[b], beta_b[b]) +
        (a0 - beta_a[b]) * u1 + (b0 - beta_b[b]) * u0;
    }
  }
  return value;
}

double nig_elbo(const NumericVector& nig_m,
                const NumericVector& nig_kappa,
                const NumericVector& nig_shape,
                const NumericVector& nig_scale,
                const double prior_mean,
                const double prior_kappa,
                const double prior_shape,
                const double prior_scale) {
  double value = 0.0;
  for (int b = 0; b < nig_m.size(); ++b) {
    const double expected_log_variance =
      std::log(nig_scale[b]) - R::digamma(nig_shape[b]);
    const double expected_precision = nig_shape[b] / nig_scale[b];
    const double mean_difference = nig_m[b] - prior_mean;

    const double expected_log_prior_variance =
      prior_shape * std::log(prior_scale) -
      R::lgammafn(prior_shape) -
      (prior_shape + 1.0) * expected_log_variance -
      prior_scale * expected_precision;
    const double expected_log_prior_mean =
      0.5 * (std::log(prior_kappa) - std::log(2.0 * M_PI) -
             expected_log_variance) -
      0.5 * prior_kappa *
        (1.0 / nig_kappa[b] +
         mean_difference * mean_difference * expected_precision);

    const double expected_log_q_variance =
      nig_shape[b] * std::log(nig_scale[b]) -
      R::lgammafn(nig_shape[b]) -
      (nig_shape[b] + 1.0) * expected_log_variance -
      nig_shape[b];
    const double expected_log_q_mean =
      0.5 * (std::log(nig_kappa[b]) - std::log(2.0 * M_PI) -
             expected_log_variance) - 0.5;

    value +=
      expected_log_prior_variance + expected_log_prior_mean -
      expected_log_q_variance - expected_log_q_mean;
  }
  return value;
}

ElboResult compute_elbo(const NumericMatrix& tau,
                        const NumericMatrix& edge_scores,
                        const int K,
                        const NumericVector& beta_a,
                        const NumericVector& beta_b,
                        const NumericVector& nig_m,
                        const NumericVector& nig_kappa,
                        const NumericVector& nig_shape,
                        const NumericVector& nig_scale,
                        const double gamma,
                        const double alpha,
                        const double beta,
                        const double alpha_prime,
                        const double beta_prime,
                        const double prior_mean,
                        const double prior_kappa,
                        const double prior_shape,
                        const double prior_scale) {
  ElboResult output;
  output.allocation = allocation_elbo(tau, gamma);
  output.edge = edge_elbo(tau, edge_scores);
  output.beta = beta_elbo(
    K, beta_a, beta_b, alpha, beta, alpha_prime, beta_prime
  );
  output.nig = nig_elbo(
    nig_m, nig_kappa, nig_shape, nig_scale,
    prior_mean, prior_kappa, prior_shape, prior_scale
  );
  return output;
}

IntegerMatrix make_edge_index(const int n) {
  IntegerMatrix index(n, n);
  std::fill(index.begin(), index.end(), -1);
  int edge = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      index(i, j) = edge;
      index(j, i) = edge;
      ++edge;
    }
  }
  return index;
}

void update_tau_sequential(NumericMatrix& tau,
                           const NumericMatrix& edge_scores,
                           const double gamma,
                           int& numerical_rebuilds,
                           const bool count_pmf_rebuild) {
  const int n = tau.nrow();
  const int K = tau.ncol();
  const IntegerMatrix edge_index = make_edge_index(n);
  std::vector<LongPmf> full_pmfs(K);
  for (int k = 0; k < K; ++k) {
    full_pmfs[k] = count_pmf_long(tau, k);
  }

  for (int i = 0; i < n; ++i) {
    // Rebuild periodically from tau rather than allowing polynomial-division
    // round-off to accumulate across an entire n-node sweep.
    if (i > 0 && i % 8 == 0) {
      for (int k = 0; k < K; ++k) {
        full_pmfs[k] = count_pmf_long(tau, k);
      }
    }

    std::vector<LongPmf> reduced_pmfs(K);
    NumericVector log_weights(K);

    for (int k = 0; k < K; ++k) {
      reduced_pmfs[k] = count_pmf_rebuild
        ? remove_factor_or_rebuild(full_pmfs[k], tau, k, i, numerical_rebuilds)
        : remove_factor_long(full_pmfs[k], tau(i, k));
      const double prior_score =
        prior_score_long(reduced_pmfs[k], gamma);

      double data_score = 0.0;
      for (int j = 0; j < n; ++j) {
        if (j == i) {
          continue;
        }
        const int edge = edge_index(i, j);
        for (int s = 0; s < K; ++s) {
          const int b = nsbmvi::block_index(k, s, K);
          data_score += tau(j, s) * edge_scores(edge, b);
        }
      }
      log_weights[k] = prior_score + data_score;
    }

    const double maximum = max(log_weights);
    double denominator = 0.0;
    for (int k = 0; k < K; ++k) {
      log_weights[k] = std::exp(log_weights[k] - maximum);
      denominator += log_weights[k];
    }
    if (!R_finite(denominator) || denominator <= 0.0) {
      stop("The tau softmax normalization failed.");
    }

    for (int k = 0; k < K; ++k) {
      tau(i, k) = log_weights[k] / denominator;
      full_pmfs[k] = add_factor_long(reduced_pmfs[k], tau(i, k));
    }
  }
}

NumericVector cumulative_qvalues(const NumericVector& lvalue) {
  const int M = lvalue.size();
  std::vector<int> order(M);
  for (int edge = 0; edge < M; ++edge) order[edge] = edge;
  std::stable_sort(order.begin(), order.end(),
    [&lvalue](const int left, const int right) {
      return lvalue[left] < lvalue[right];
    });

  NumericVector qvalue(M);
  double cumulative = 0.0;
  for (int rank = 0; rank < M; ++rank) {
    cumulative += lvalue[order[rank]];
    qvalue[order[rank]] = cumulative / static_cast<double>(rank + 1);
  }
  return qvalue;
}

List averaged_graph_output(const NumericMatrix& tau,
                           const NumericMatrix& phi) {
  const NumericMatrix weights = nsbmvi::compute_block_weights(tau);
  const int M = weights.nrow();
  NumericVector lvalue(M);

  for (int edge = 0; edge < M; ++edge) {
    for (int b = 0; b < weights.ncol(); ++b) {
      lvalue[edge] += weights(edge, b) * (1.0 - phi(edge, b));
    }
  }

  return List::create(
    Named("lvalue") = lvalue,
    Named("qvalue") = cumulative_qvalues(lvalue)
  );
}

std::pair<double, double> gaussian_rejection_probabilities(
    const double threshold,
    const double edge_probability,
    const double signal_mean,
    const double signal_sd) {
  if (threshold >= 1.0) return std::make_pair(1.0, 1.0);
  if (!(threshold > 0.0) || !(edge_probability > 0.0) ||
      !(edge_probability < 1.0) || !(signal_sd > 0.0)) {
    return std::make_pair(0.0, 0.0);
  }

  const double a = 1.0 / (signal_sd * signal_sd) - 1.0;
  const double b = -2.0 * signal_mean / (signal_sd * signal_sd);
  const double c = signal_mean * signal_mean /
    (signal_sd * signal_sd) + 2.0 * std::log(
      signal_sd * (1.0 / edge_probability - 1.0) *
      (1.0 / threshold - 1.0)
    );
  double null_probability = 0.0;
  double signal_probability = 0.0;

  if (std::abs(a) > 1e-12) {
    if (a < 0.0) {
      null_probability = 1.0;
      signal_probability = 1.0;
    }
    const double discriminant = b * b - 4.0 * a * c;
    if (discriminant > 0.0) {
      const double root = std::sqrt(discriminant);
      const double first = (-b + root) / (2.0 * a);
      const double second = (-b - root) / (2.0 * a);
      null_probability +=
        R::pnorm(first, 0.0, 1.0, true, false) -
        R::pnorm(second, 0.0, 1.0, true, false);
      signal_probability +=
        R::pnorm(first, signal_mean, signal_sd, true, false) -
        R::pnorm(second, signal_mean, signal_sd, true, false);
    }
  } else if (std::abs(b) > 1e-12) {
    const double cutoff = -c / b;
    if (b < 0.0) {
      null_probability = 1.0 -
        R::pnorm(cutoff, 0.0, 1.0, true, false);
      signal_probability = 1.0 -
        R::pnorm(cutoff, signal_mean, signal_sd, true, false);
    } else {
      null_probability = R::pnorm(cutoff, 0.0, 1.0, true, false);
      signal_probability = R::pnorm(
        cutoff, signal_mean, signal_sd, true, false
      );
    }
  } else {
    const double value = threshold >= 1.0 - edge_probability ? 1.0 : 0.0;
    null_probability = value;
    signal_probability = value;
  }

  null_probability = std::min(std::max(null_probability, 0.0), 1.0);
  signal_probability = std::min(std::max(signal_probability, 0.0), 1.0);
  return std::make_pair(null_probability, signal_probability);
}

List plugin_graph_output(const NumericMatrix& X,
                         const NumericMatrix& tau,
                         const NumericVector& beta_a,
                         const NumericVector& beta_b,
                         const NumericVector& nig_m,
                         const NumericVector& nig_shape,
                         const NumericVector& nig_scale) {
  const int n = X.nrow();
  const int K = tau.ncol();
  const int B = nsbmvi::block_count(K);
  const int M = n * (n - 1) / 2;
  IntegerVector labels(n);
  NumericVector proportions(K);
  for (int i = 0; i < n; ++i) {
    int best = 0;
    for (int k = 1; k < K; ++k) {
      if (tau(i, k) > tau(i, best)) best = k;
    }
    labels[i] = best + 1;
    proportions[best] += 1.0 / n;
  }

  NumericVector edge_probability(B);
  NumericVector signal_mean(B);
  NumericVector signal_variance(B);
  NumericVector signal_sd(B);
  for (int b = 0; b < B; ++b) {
    edge_probability[b] = beta_a[b] / (beta_a[b] + beta_b[b]);
    signal_mean[b] = nig_m[b];
    if (nig_shape[b] <= 1.0) {
      stop("The plug-in posterior mean of sigma2 requires nig_shape > 1.");
    }
    signal_variance[b] = nig_scale[b] / (nig_shape[b] - 1.0);
    signal_sd[b] = std::sqrt(signal_variance[b]);
  }

  NumericVector lvalue(M);
  int edge = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const int b = nsbmvi::block_index(labels[i] - 1, labels[j] - 1, K);
      const double q = std::min(std::max(edge_probability[b], 1e-12),
                                1.0 - 1e-12);
      const double log_null = std::log1p(-q) +
        R::dnorm(X(i, j), 0.0, 1.0, true);
      const double log_signal = std::log(q) +
        R::dnorm(X(i, j), signal_mean[b], signal_sd[b], true);
      lvalue[edge] = nsbmvi::stable_logistic(log_null - log_signal);
      ++edge;
    }
  }

  NumericVector mfdr_numerator(M);
  NumericVector mfdr_denominator(M);
  for (int r = 0; r < K; ++r) {
    for (int s = r; s < K; ++s) {
      const int b = nsbmvi::block_index(r, s, K);
      const double symmetry = r == s ? 1.0 : 2.0;
      const double block_weight = symmetry * proportions[r] * proportions[s];
      const double null_weight = block_weight * (1.0 - edge_probability[b]);
      const double signal_weight = block_weight * edge_probability[b];
      for (int current = 0; current < M; ++current) {
        const std::pair<double, double> probabilities =
          gaussian_rejection_probabilities(
            lvalue[current], edge_probability[b],
            signal_mean[b], signal_sd[b]
          );
        mfdr_numerator[current] += null_weight * probabilities.first;
        mfdr_denominator[current] += signal_weight * probabilities.second;
      }
    }
  }

  NumericVector mfdr_qvalue(M);
  for (int current = 0; current < M; ++current) {
    const double total = mfdr_numerator[current] + mfdr_denominator[current];
    mfdr_qvalue[current] = total > 0.0
      ? mfdr_numerator[current] / total
      : 0.0;
  }

  return List::create(
    Named("lvalue") = lvalue,
    Named("qvalue_cumulative") = cumulative_qvalues(lvalue),
    Named("qvalue_mfdr") = mfdr_qvalue,
    Named("labels") = labels,
    Named("pi") = proportions,
    Named("Q") = edge_probability,
    Named("mu") = signal_mean,
    Named("sigma2") = signal_variance
  );
}

}  // namespace


// [[Rcpp::export]]
List nsbm_vi_fixed_k_cpp(
    NumericMatrix X,
    const int K,
    Nullable<NumericMatrix> tau_init = R_NilValue,
    Nullable<List> hyperparameters = R_NilValue,
    const int max_iter = 300,
    const double tolerance = 1e-6,
    const int min_iter = 2,
    const bool verbose = false,
    Nullable<NumericMatrix> phi_init = R_NilValue,
    const bool count_pmf_rebuild = true) {
  nsbmvi::validate_X(X);
  if (K < 1) {
    stop("K must be positive.");
  }
  if (max_iter < 1 || min_iter < 1 || min_iter > max_iter) {
    stop("Require 1 <= min_iter <= max_iter.");
  }
  if (!R_finite(tolerance) || tolerance <= 0.0) {
    stop("tolerance must be positive and finite.");
  }

  const int n = X.nrow();
  const int B = nsbmvi::block_count(K);
  NumericMatrix tau;
  if (tau_init.isNotNull()) {
    tau = as<NumericMatrix>(tau_init);
    if (tau.nrow() != n || tau.ncol() != K) {
      stop("tau_init must be an n by K matrix.");
    }
    tau = clone(tau);
  } else {
    tau = initialize_tau(n, K);
  }
  nsbmvi::validate_tau(tau);

  List hp = hyperparameters.isNotNull()
    ? as<List>(hyperparameters)
    : List::create();
  const double alpha = get_hyperparameter(hp, "alpha", 2.0);
  const double beta = get_hyperparameter(hp, "beta", 2.0);
  const double alpha_prime = get_hyperparameter(hp, "alpha_prime", 2.0);
  const double beta_prime = get_hyperparameter(hp, "beta_prime", 2.0);
  const double prior_mean = get_hyperparameter(hp, "rou", 0.0);
  const double prior_kappa = get_hyperparameter(hp, "kappa", 1.0);
  const double prior_shape = get_hyperparameter(hp, "delta", 3.0);
  const double prior_scale = get_hyperparameter(hp, "xi", 2.0);
  const double gamma = get_hyperparameter(hp, "gamma", 1.0);

  if (alpha <= 0.0 || beta <= 0.0 ||
      alpha_prime <= 0.0 || beta_prime <= 0.0 ||
      prior_kappa <= 0.0 || prior_shape <= 0.0 ||
      prior_scale <= 0.0 || gamma <= 0.0) {
    stop("All Beta, NIG, and gamma hyperparameters must be positive.");
  }

  NumericVector beta_a(B);
  NumericVector beta_b(B);
  NumericVector nig_m(B, prior_mean);
  NumericVector nig_kappa(B, prior_kappa);
  NumericVector nig_shape(B, prior_shape);
  NumericVector nig_scale(B, prior_scale);
  for (int r = 0; r < K; ++r) {
    for (int s = r; s < K; ++s) {
      const int b = nsbmvi::block_index(r, s, K);
      beta_a[b] = (r == s) ? alpha_prime : alpha;
      beta_b[b] = (r == s) ? beta_prime : beta;
    }
  }

  NumericMatrix phi;
  if (phi_init.isNotNull()) {
    phi = as<NumericMatrix>(phi_init);
    const int M = n * (n - 1) / 2;
    if (phi.nrow() != M || phi.ncol() != B) {
      stop("phi_init must have n * (n - 1) / 2 rows and K * (K + 1) / 2 columns.");
    }
    phi = clone(phi);
    for (int edge = 0; edge < M; ++edge) {
      for (int b = 0; b < B; ++b) {
        if (!R_finite(phi(edge, b)) ||
            phi(edge, b) < 0.0 || phi(edge, b) > 1.0) {
          stop("Every phi_init entry must be finite and lie in [0, 1].");
        }
      }
    }
  } else {
    phi = nsbmvi::update_phi(
      X, K, beta_a, beta_b,
      nig_m, nig_kappa, nig_shape, nig_scale
    );
  }
  List posterior = nsbmvi::update_block_posteriors(
    X, tau, phi,
    alpha, beta, alpha_prime, beta_prime,
    prior_mean, prior_kappa, prior_shape, prior_scale
  );
  beta_a = posterior["beta_a"];
  beta_b = posterior["beta_b"];
  nig_m = posterior["nig_m"];
  nig_kappa = posterior["nig_kappa"];
  nig_shape = posterior["nig_shape"];
  nig_scale = posterior["nig_scale"];

  NumericMatrix edge_scores = nsbmvi::compute_edge_scores(
    X, phi, K, beta_a, beta_b,
    nig_m, nig_kappa, nig_shape, nig_scale
  );
  ElboResult current = compute_elbo(
    tau, edge_scores, K,
    beta_a, beta_b, nig_m, nig_kappa, nig_shape, nig_scale,
    gamma, alpha, beta, alpha_prime, beta_prime,
    prior_mean, prior_kappa, prior_shape, prior_scale
  );

  NumericVector elbo_trace(max_iter + 1);
  elbo_trace[0] = current.total();
  NumericMatrix component_trace(max_iter + 1, 4);
  component_trace(0, 0) = current.allocation;
  component_trace(0, 1) = current.edge;
  component_trace(0, 2) = current.beta;
  component_trace(0, 3) = current.nig;
  colnames(component_trace) = CharacterVector::create(
    "allocation", "edge", "beta", "nig"
  );

  bool converged = false;
  bool monotone = true;
  int completed_iterations = 0;
  int numerical_rebuilds = 0;

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    phi = nsbmvi::update_phi(
      X, K, beta_a, beta_b,
      nig_m, nig_kappa, nig_shape, nig_scale
    );
    posterior = nsbmvi::update_block_posteriors(
      X, tau, phi,
      alpha, beta, alpha_prime, beta_prime,
      prior_mean, prior_kappa, prior_shape, prior_scale
    );
    beta_a = posterior["beta_a"];
    beta_b = posterior["beta_b"];
    nig_m = posterior["nig_m"];
    nig_kappa = posterior["nig_kappa"];
    nig_shape = posterior["nig_shape"];
    nig_scale = posterior["nig_scale"];

    edge_scores = nsbmvi::compute_edge_scores(
      X, phi, K, beta_a, beta_b,
      nig_m, nig_kappa, nig_shape, nig_scale
    );
    update_tau_sequential(tau, edge_scores, gamma, numerical_rebuilds, count_pmf_rebuild);

    current = compute_elbo(
      tau, edge_scores, K,
      beta_a, beta_b, nig_m, nig_kappa, nig_shape, nig_scale,
      gamma, alpha, beta, alpha_prime, beta_prime,
      prior_mean, prior_kappa, prior_shape, prior_scale
    );
    elbo_trace[iteration] = current.total();
    component_trace(iteration, 0) = current.allocation;
    component_trace(iteration, 1) = current.edge;
    component_trace(iteration, 2) = current.beta;
    component_trace(iteration, 3) = current.nig;
    completed_iterations = iteration;

    const double improvement =
      elbo_trace[iteration] - elbo_trace[iteration - 1];
    const double relative_improvement =
      std::abs(improvement) /
      (1.0 + std::abs(elbo_trace[iteration - 1]));
    if (improvement < -1e-7 *
        (1.0 + std::abs(elbo_trace[iteration - 1]))) {
      monotone = false;
    }

    if (verbose) {
      Rcout << "iteration=" << iteration
            << " elbo=" << elbo_trace[iteration]
            << " change=" << improvement << "\n";
    }

    if (iteration >= min_iter && relative_improvement < tolerance) {
      converged = true;
      break;
    }
  }

  elbo_trace = elbo_trace[Range(0, completed_iterations)];
  component_trace = component_trace(
    Range(0, completed_iterations), Range(0, 3)
  );
  // Refresh q(A | Z) using the final conjugate block factors.  During a CAVI
  // sweep phi is updated before those factors, so the stored phi otherwise
  // lags the final parameter update by one coordinate step.
  phi = nsbmvi::update_phi(
    X, K, beta_a, beta_b,
    nig_m, nig_kappa, nig_shape, nig_scale
  );
  const List averaged_graph = averaged_graph_output(tau, phi);
  const List plugin_graph = plugin_graph_output(
    X, tau, beta_a, beta_b, nig_m, nig_shape, nig_scale
  );

  return List::create(
    Named("K") = K,
    Named("tau") = tau,
    Named("phi") = phi,
    Named("beta_a") = beta_a,
    Named("beta_b") = beta_b,
    Named("nig_m") = nig_m,
    Named("nig_kappa") = nig_kappa,
    Named("nig_shape") = nig_shape,
    Named("nig_scale") = nig_scale,
    Named("elbo_trace") = elbo_trace,
    Named("elbo_components") = component_trace,
    // Backward-compatible aliases: lvalue/qvalue are the uncertainty-averaged
    // variational marginals and their cumulative posterior-FDP q-values.
    Named("lvalue") = averaged_graph["lvalue"],
    Named("qvalue") = averaged_graph["qvalue"],
    Named("averaged_lvalue") = averaged_graph["lvalue"],
    Named("averaged_qvalue") = averaged_graph["qvalue"],
    Named("plugin_lvalue") = plugin_graph["lvalue"],
    Named("plugin_qvalue_cumulative") = plugin_graph["qvalue_cumulative"],
    Named("plugin_qvalue_mfdr") = plugin_graph["qvalue_mfdr"],
    Named("plugin_parameters") = List::create(
      Named("labels") = plugin_graph["labels"],
      Named("pi") = plugin_graph["pi"],
      Named("Q") = plugin_graph["Q"],
      Named("mu") = plugin_graph["mu"],
      Named("sigma2") = plugin_graph["sigma2"]
    ),
    Named("converged") = converged,
    Named("monotone") = monotone,
    Named("iterations") = completed_iterations,
    Named("count_pmf_rebuilds") = numerical_rebuilds,
    Named("hyperparameters") = List::create(
      Named("alpha") = alpha,
      Named("beta") = beta,
      Named("alpha_prime") = alpha_prime,
      Named("beta_prime") = beta_prime,
      Named("rou") = prior_mean,
      Named("kappa") = prior_kappa,
      Named("delta") = prior_shape,
      Named("xi") = prior_scale,
      Named("gamma") = gamma
    )
  );
}
