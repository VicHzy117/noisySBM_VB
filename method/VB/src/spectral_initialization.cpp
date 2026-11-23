// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double LOG_TWO_PI = 1.83787706640934548356;

struct EmFit {
  bool valid = false;
  bool converged = false;
  bool signal_detected = false;
  bool two_sided = false;
  double q = 0.0;
  double mu = 0.0;
  double sigma2 = 1.0;
  double eff_n = 0.0;
  double loglik = NA_REAL;
  double loglik_null = NA_REAL;
  double lr_stat = NA_REAL;
  int iterations = 0;
};

double clamp_probability(const double value, const double floor) {
  return std::min(std::max(value, floor), 1.0 - floor);
}

double log_normal_density(const double x, const double mean,
                          const double variance) {
  return -0.5 * (LOG_TWO_PI + std::log(variance) +
                 (x - mean) * (x - mean) / variance);
}

double logistic(const double value) {
  if (value >= 0.0) {
    return 1.0 / (1.0 + std::exp(-value));
  }
  const double exponential = std::exp(value);
  return exponential / (1.0 + exponential);
}

double log_add_exp(const double left, const double right) {
  const double maximum = std::max(left, right);
  return maximum + std::log(std::exp(left - maximum) +
                            std::exp(right - maximum));
}

double mixture_loglik(const std::vector<double>& x, const double q,
                      const double mu, const double sigma2,
                      const double null_variance) {
  const double log_q = std::log(q);
  const double log_one_minus_q = std::log1p(-q);
  double value = 0.0;
  for (const double observation : x) {
    const double signal = log_q +
      log_normal_density(observation, mu, sigma2);
    const double null_value = log_one_minus_q +
      log_normal_density(observation, 0.0, null_variance);
    value += log_add_exp(signal, null_value);
  }
  return value;
}

double null_loglik(const std::vector<double>& x,
                   const double null_variance) {
  double value = 0.0;
  for (const double observation : x) {
    value += log_normal_density(observation, 0.0, null_variance);
  }
  return value;
}

double type7_quantile(const std::vector<double>& sorted,
                      const double probability) {
  if (sorted.size() == 1U) return sorted.front();
  const double index = (sorted.size() - 1.0) * probability;
  const std::size_t lower = static_cast<std::size_t>(std::floor(index));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(index));
  const double fraction = index - lower;
  return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper];
}

EmFit run_em_start(const std::vector<double>& x, const double mu_start,
                   const double null_variance, const int max_iter,
                   const double tolerance, const double min_eff_n,
                   const double sigma2_floor,
                   const double probability_floor) {
  EmFit fit;
  double q = 0.2;
  double mu = mu_start;
  double sigma2 = std::max(null_variance, sigma2_floor);
  double old_loglik = -std::numeric_limits<double>::infinity();
  double effective_n = 0.0;

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    const double log_q = std::log(q);
    const double log_one_minus_q = std::log1p(-q);
    double loglik = 0.0;
    double weighted_sum = 0.0;
    double weighted_square_sum = 0.0;
    effective_n = 0.0;

    for (const double observation : x) {
      const double signal = log_q +
        log_normal_density(observation, mu, sigma2);
      const double null_value = log_one_minus_q +
        log_normal_density(observation, 0.0, null_variance);
      const double probability = logistic(signal - null_value);
      loglik += log_add_exp(signal, null_value);
      effective_n += probability;
      weighted_sum += probability * observation;
      weighted_square_sum += probability * observation * observation;
    }

    if (!std::isfinite(effective_n) || effective_n < min_eff_n) {
      return fit;
    }

    q = clamp_probability(effective_n / x.size(), probability_floor);
    mu = weighted_sum / effective_n;
    sigma2 = std::max(
      weighted_square_sum / effective_n - mu * mu,
      sigma2_floor
    );

    fit.iterations = iteration;
    if (std::isfinite(loglik) &&
        std::abs(loglik - old_loglik) < tolerance) {
      fit.converged = true;
      break;
    }
    old_loglik = loglik;
  }

  fit.valid = true;
  fit.q = q;
  fit.mu = mu;
  fit.sigma2 = sigma2;
  fit.eff_n = effective_n;
  fit.loglik = mixture_loglik(
    x, q, mu, sigma2, null_variance
  );
  return fit;
}

EmFit pooled_em(const std::vector<double>& x, const double null_sd,
                const int n_starts, const int max_iter,
                const double tolerance, const double min_eff_n,
                const double sigma2_floor,
                const double probability_floor) {
  std::vector<double> sorted = x;
  std::sort(sorted.begin(), sorted.end());
  const double null_variance = null_sd * null_sd;
  EmFit best;

  for (int start = 0; start < n_starts; ++start) {
    const double probability = n_starts == 1
      ? 0.60
      : 0.60 + (0.99 - 0.60) * start / (n_starts - 1.0);
    EmFit candidate = run_em_start(
      x, type7_quantile(sorted, probability), null_variance,
      max_iter, tolerance, min_eff_n, sigma2_floor, probability_floor
    );
    if (candidate.valid &&
        (!best.valid || candidate.loglik > best.loglik)) {
      best = candidate;
    }
  }

  if (!best.valid) {
    best.q = 0.0;
    best.mu = 0.0;
    best.sigma2 = null_variance;
    return best;
  }

  best.loglik_null = null_loglik(x, null_variance);
  best.lr_stat = 2.0 * (best.loglik - best.loglik_null);
  if (!std::isfinite(best.lr_stat) ||
      best.lr_stat < 3.0 * std::log(static_cast<double>(x.size()))) {
    best.signal_detected = false;
    return best;
  }

  best.signal_detected = true;
  const double cutoff = 1.96 * null_sd;
  int lower_count = 0;
  int upper_count = 0;
  for (const double observation : x) {
    lower_count += observation < -cutoff;
    upper_count += observation > cutoff;
  }
  const double expected_tail = x.size() * (1.0 - best.q) *
    R::pnorm(-1.96, 0.0, 1.0, true, false);
  best.two_sided =
    lower_count > 20 && lower_count > 1.5 * expected_tail &&
    upper_count > 20 && upper_count > 1.5 * expected_tail;
  return best;
}

std::vector<double> upper_triangle_rowwise(const NumericMatrix& X) {
  const int n = X.nrow();
  std::vector<double> values;
  values.reserve(n * (n - 1) / 2);
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      values.push_back(X(i, j));
    }
  }
  return values;
}

arma::mat pooled_probability_matrix(const NumericMatrix& X,
                                    const EmFit& fit,
                                    const double null_variance,
                                    const double probability_floor) {
  const int n = X.nrow();
  arma::mat probabilities(n, n, arma::fill::zeros);
  const double q = clamp_probability(fit.q, probability_floor);
  const double log_q = std::log(q);
  const double log_one_minus_q = std::log1p(-q);
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double observation = X(i, j);
      const double signal = log_q +
        log_normal_density(observation, fit.mu, fit.sigma2);
      const double null_value = log_one_minus_q +
        log_normal_density(observation, 0.0, null_variance);
      const double probability = clamp_probability(
        logistic(signal - null_value), probability_floor
      );
      probabilities(i, j) = probability;
      probabilities(j, i) = probability;
    }
  }
  return probabilities;
}

double squared_distance(const arma::rowvec& left,
                        const arma::rowvec& right) {
  return arma::accu(arma::square(left - right));
}

std::vector<int> kmeans_once(const arma::mat& embedding, const int K,
                             const int max_iter,
                             std::mt19937& generator,
                             double& within_sum_squares) {
  const int n = embedding.n_rows;
  const int dimension = embedding.n_cols;
  arma::mat centers(K, dimension, arma::fill::zeros);
  std::vector<int> center_indices;
  center_indices.reserve(K);
  std::uniform_int_distribution<int> first_distribution(0, n - 1);
  center_indices.push_back(first_distribution(generator));
  centers.row(0) = embedding.row(center_indices[0]);

  std::vector<double> minimum_distance(n,
    std::numeric_limits<double>::infinity());
  for (int center = 1; center < K; ++center) {
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
      minimum_distance[i] = std::min(
        minimum_distance[i],
        squared_distance(embedding.row(i), centers.row(center - 1))
      );
      total += minimum_distance[i];
    }

    int selected = -1;
    if (total > 0.0 && std::isfinite(total)) {
      std::uniform_real_distribution<double> draw(0.0, total);
      const double target = draw(generator);
      double cumulative = 0.0;
      for (int i = 0; i < n; ++i) {
        cumulative += minimum_distance[i];
        if (cumulative >= target) {
          selected = i;
          break;
        }
      }
    }
    if (selected < 0 ||
        std::find(center_indices.begin(), center_indices.end(), selected) !=
          center_indices.end()) {
      for (int i = 0; i < n; ++i) {
        if (std::find(center_indices.begin(), center_indices.end(), i) ==
            center_indices.end()) {
          selected = i;
          break;
        }
      }
    }
    center_indices.push_back(selected);
    centers.row(center) = embedding.row(selected);
  }

  std::vector<int> labels(n, -1);
  std::vector<int> previous(n, -2);
  std::vector<double> assigned_distance(n, 0.0);

  for (int iteration = 0; iteration < max_iter; ++iteration) {
    previous = labels;
    std::vector<int> counts(K, 0);
    for (int i = 0; i < n; ++i) {
      double best_distance = std::numeric_limits<double>::infinity();
      int best_cluster = 0;
      for (int cluster = 0; cluster < K; ++cluster) {
        const double distance = squared_distance(
          embedding.row(i), centers.row(cluster)
        );
        if (distance < best_distance) {
          best_distance = distance;
          best_cluster = cluster;
        }
      }
      labels[i] = best_cluster;
      assigned_distance[i] = best_distance;
      ++counts[best_cluster];
    }

    for (int cluster = 0; cluster < K; ++cluster) {
      if (counts[cluster] != 0) continue;
      int replacement = -1;
      double largest_distance = -1.0;
      for (int i = 0; i < n; ++i) {
        if (counts[labels[i]] > 1 && assigned_distance[i] > largest_distance) {
          replacement = i;
          largest_distance = assigned_distance[i];
        }
      }
      if (replacement >= 0) {
        --counts[labels[replacement]];
        labels[replacement] = cluster;
        ++counts[cluster];
      }
    }

    centers.zeros();
    for (int i = 0; i < n; ++i) {
      centers.row(labels[i]) += embedding.row(i);
    }
    for (int cluster = 0; cluster < K; ++cluster) {
      centers.row(cluster) /= counts[cluster];
    }
    if (labels == previous) break;
  }

  within_sum_squares = 0.0;
  for (int i = 0; i < n; ++i) {
    within_sum_squares += squared_distance(
      embedding.row(i), centers.row(labels[i])
    );
  }
  return labels;
}

std::vector<int> multistart_kmeans(const arma::mat& embedding, const int K,
                                   const int n_starts, const int seed) {
  if (K == 1) return std::vector<int>(embedding.n_rows, 0);
  std::mt19937 generator(static_cast<unsigned int>(seed));
  std::vector<int> best;
  double best_within = std::numeric_limits<double>::infinity();
  for (int start = 0; start < n_starts; ++start) {
    double within = 0.0;
    std::vector<int> candidate = kmeans_once(
      embedding, K, 100, generator, within
    );
    if (within < best_within) {
      best_within = within;
      best = candidate;
    }
  }
  return best;
}

NumericMatrix soften_labels(const std::vector<int>& labels, const int K,
                            const double epsilon) {
  const int n = labels.size();
  NumericMatrix tau(n, K);
  if (K == 1) {
    std::fill(tau.begin(), tau.end(), 1.0);
    return tau;
  }
  const double other = epsilon / (K - 1.0);
  std::fill(tau.begin(), tau.end(), other);
  for (int i = 0; i < n; ++i) tau(i, labels[i]) = 1.0 - epsilon;
  return tau;
}

bool all_groups_occupied(const std::vector<int>& labels, const int K) {
  std::vector<int> counts(K, 0);
  for (const int label : labels) ++counts[label];
  return std::all_of(counts.begin(), counts.end(),
                     [](const int count) { return count > 0; });
}

std::vector<int> perturb_labels(const std::vector<int>& labels, const int K,
                                const double fraction,
                                std::mt19937& generator) {
  if (K == 1 || fraction <= 0.0) return labels;
  const int n = labels.size();
  const int n_change = std::max(1, static_cast<int>(std::round(n * fraction)));
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), generator);
  std::uniform_int_distribution<int> alternative(0, K - 2);
  std::vector<int> output = labels;
  for (int position = 0; position < n_change; ++position) {
    const int i = indices[position];
    int new_label = alternative(generator);
    if (new_label >= labels[i]) ++new_label;
    output[i] = new_label;
  }
  return all_groups_occupied(output, K) ? output : labels;
}

IntegerVector one_based_labels(const std::vector<int>& labels) {
  IntegerVector output(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i) output[i] = labels[i] + 1;
  return output;
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();
}

}  // namespace


// [[Rcpp::export]]
List nsbm_spectral_init_cpp(
    NumericMatrix X,
    const int K,
    const int seed = 1,
    const double epsilon = 0.1,
    const int kmeans_nstart = 25,
    const int n_perturb = 2,
    const double perturb_fraction = 0.15,
    const double null_sd = 1.0,
    const int em_n_starts = 5,
    const int em_max_iter = 5000,
    const double em_tolerance = 1e-8,
    const double min_eff_n = 50.0,
    const double sigma2_floor = 1e-3,
    const double probability_floor = 1e-8) {
  const int n = X.nrow();
  if (n != X.ncol() || n < 2) stop("X must be a square matrix with n >= 2.");
  if (K < 1 || K >= n) stop("K must be between 1 and n - 1.");
  if (epsilon <= 0.0 || epsilon >= 1.0) stop("epsilon must lie in (0, 1).");
  if (kmeans_nstart < 1 || n_perturb < 0 ||
      perturb_fraction < 0.0 || perturb_fraction > 1.0) {
    stop("Invalid k-means or perturbation control.");
  }
  if (null_sd <= 0.0 || em_n_starts < 1 || em_max_iter < 1 ||
      em_tolerance <= 0.0 || min_eff_n <= 0.0 || sigma2_floor <= 0.0 ||
      probability_floor <= 0.0 || probability_floor >= 0.5) {
    stop("Invalid pooled-EM control.");
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (!R_finite(X(i, j))) stop("X must contain only finite values.");
      if (std::abs(X(i, j) - X(j, i)) > 1e-8) stop("X must be symmetric.");
    }
  }

  const std::vector<double> observations = upper_triangle_rowwise(X);
  const auto pooled_start = std::chrono::steady_clock::now();
  EmFit pooled = pooled_em(
    observations, null_sd, em_n_starts, em_max_iter, em_tolerance,
    min_eff_n, sigma2_floor, probability_floor
  );
  const double pooled_time = elapsed_seconds(pooled_start);

  const auto spectral_start = std::chrono::steady_clock::now();
  arma::mat edge_probability(n, n, arma::fill::zeros);
  arma::mat spectral_weights(n, n, arma::fill::zeros);
  std::string method;
  if (pooled.signal_detected && !pooled.two_sided) {
    edge_probability = pooled_probability_matrix(
      X, pooled, null_sd * null_sd, probability_floor
    );
    spectral_weights = edge_probability;
    method = "pooled-posterior spectral (Rcpp)";
  } else {
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double value = X(i, j) * X(i, j) - 1.0;
        spectral_weights(i, j) = value;
        spectral_weights(j, i) = value;
        edge_probability(i, j) = 0.1;
        edge_probability(j, i) = 0.1;
      }
    }
    method = "squared-residual fallback (Rcpp)";
  }

  const double off_diagonal_mean = arma::accu(arma::trimatu(spectral_weights, 1)) /
    (n * (n - 1.0) / 2.0);
  arma::mat centered = spectral_weights;
  for (int i = 0; i < n; ++i) {
    centered(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      centered(i, j) -= off_diagonal_mean;
      centered(j, i) = centered(i, j);
    }
  }

  arma::vec eigenvalues_all;
  arma::mat eigenvectors_all;
  if (!arma::eig_sym(eigenvalues_all, eigenvectors_all, centered)) {
    stop("The symmetric eigendecomposition failed.");
  }
  std::vector<std::pair<double, int>> eigen_order(n);
  for (int index = 0; index < n; ++index) {
    eigen_order[index] = std::make_pair(
      std::abs(eigenvalues_all[index]), index
    );
  }
  std::sort(eigen_order.begin(), eigen_order.end(),
            [](const std::pair<double, int>& left,
               const std::pair<double, int>& right) {
              return left.first > right.first;
            });

  arma::mat embedding(n, K);
  NumericVector selected_eigenvalues(K);
  IntegerVector selected_indices(K);
  for (int component = 0; component < K; ++component) {
    const int index = eigen_order[component].second;
    embedding.col(component) = eigenvectors_all.col(index);
    selected_eigenvalues[component] = eigenvalues_all[index];
    selected_indices[component] = index + 1;
  }
  for (int i = 0; i < n; ++i) {
    const double norm = arma::norm(embedding.row(i), 2);
    if (norm > std::sqrt(std::numeric_limits<double>::epsilon())) {
      embedding.row(i) /= norm;
    }
  }

  const std::vector<int> labels = multistart_kmeans(
    embedding, K, kmeans_nstart, seed
  );
  std::mt19937 perturb_generator(static_cast<unsigned int>(seed + 7919));
  List starts(n_perturb + 1);
  CharacterVector start_names(n_perturb + 1);
  starts[0] = List::create(
    Named("labels") = one_based_labels(labels),
    Named("tau") = soften_labels(labels, K, epsilon)
  );
  start_names[0] = "spectral";
  for (int index = 0; index < n_perturb; ++index) {
    const std::vector<int> perturbed = perturb_labels(
      labels, K, perturb_fraction, perturb_generator
    );
    starts[index + 1] = List::create(
      Named("labels") = one_based_labels(perturbed),
      Named("tau") = soften_labels(perturbed, K, epsilon)
    );
    start_names[index + 1] = "spectral_perturb_" +
      std::to_string(index + 1);
  }
  starts.attr("names") = start_names;

  const int M = observations.size();
  const int B = K * (K + 1) / 2;
  NumericMatrix phi(M, B);
  int edge = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      for (int block = 0; block < B; ++block) {
        phi(edge, block) = edge_probability(i, j);
      }
      ++edge;
    }
  }
  const double spectral_time = elapsed_seconds(spectral_start);

  return List::create(
    Named("method") = method,
    Named("starts") = starts,
    Named("phi") = phi,
    Named("edge_probability") = edge_probability,
    Named("pooled_fit") = List::create(
      Named("mu") = pooled.mu,
      Named("sigma2") = pooled.sigma2,
      Named("Q") = pooled.q,
      Named("eff_n") = pooled.eff_n,
      Named("loglik") = pooled.loglik,
      Named("loglik_null") = pooled.loglik_null,
      Named("lr_stat") = pooled.lr_stat,
      Named("iter") = pooled.iterations,
      Named("converged") = pooled.converged,
      Named("signal_detected") = pooled.signal_detected,
      Named("two_sided") = pooled.two_sided
    ),
    Named("eigenvalues") = selected_eigenvalues,
    Named("selected_eigen_indices") = selected_indices,
    Named("timing") = NumericVector::create(
      Named("pooled") = pooled_time,
      Named("spectral") = spectral_time,
      Named("total") = pooled_time + spectral_time
    )
  );
}
