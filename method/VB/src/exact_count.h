#ifndef NSBM_VI_EXACT_COUNT_H
#define NSBM_VI_EXACT_COUNT_H

#include <Rcpp.h>

namespace nsbmvi {

// Polynomial division is exact algebraically, but coefficients in the
// far tails of an n around 200 Poisson-binomial pmf can accumulate absolute
// round-off around 1e-10. Values within this tolerance are clipped before
// renormalization; larger negatives still indicate a genuine numerical fault.
constexpr double kProbabilityTolerance = 1e-8;

inline double checked_probability(const double value, const char* name) {
  if (!R_finite(value) || value < 0.0 || value > 1.0) {
    Rcpp::stop("%s must be finite and lie in [0, 1].", name);
  }
  return value;
}

inline Rcpp::NumericVector normalize_pmf(Rcpp::NumericVector pmf) {
  double total = 0.0;

  for (R_xlen_t i = 0; i < pmf.size(); ++i) {
    if (!R_finite(pmf[i])) {
      Rcpp::stop("Count pmf contains a non-finite value.");
    }
    if (pmf[i] < -kProbabilityTolerance) {
      Rcpp::stop(
        "Count pmf contains a materially negative value: %.17g at index %d.",
        pmf[i], static_cast<int>(i)
      );
    }
    if (pmf[i] < 0.0) {
      pmf[i] = 0.0;
    }
    total += pmf[i];
  }

  if (!R_finite(total) || total <= 0.0) {
    Rcpp::stop("Count pmf has non-positive or non-finite total mass.");
  }

  for (R_xlen_t i = 0; i < pmf.size(); ++i) {
    pmf[i] /= total;
  }

  return pmf;
}

inline Rcpp::NumericVector count_pmf(
    const Rcpp::NumericVector& probabilities) {
  const int n = probabilities.size();
  Rcpp::NumericVector pmf(n + 1);
  pmf[0] = 1.0;

  int current_degree = 0;
  for (int i = 0; i < n; ++i) {
    const double p = checked_probability(probabilities[i], "probability");
    Rcpp::NumericVector next(n + 1);

    for (int count = 0; count <= current_degree; ++count) {
      next[count] += (1.0 - p) * pmf[count];
      next[count + 1] += p * pmf[count];
    }

    pmf = next;
    ++current_degree;
  }

  return normalize_pmf(pmf);
}

inline Rcpp::NumericVector remove_bernoulli_factor(
    const Rcpp::NumericVector& full_pmf,
    const double probability) {
  const double p = checked_probability(probability, "probability");
  const int full_size = full_pmf.size();

  if (full_size < 2) {
    Rcpp::stop("full_pmf must have length at least two.");
  }

  Rcpp::NumericVector full = Rcpp::clone(full_pmf);
  full = normalize_pmf(full);

  const int reduced_size = full_size - 1;
  Rcpp::NumericVector reduced(reduced_size);

  if (p <= 0.5) {
    const double denominator = 1.0 - p;
    reduced[0] = full[0] / denominator;
    for (int count = 1; count < reduced_size; ++count) {
      reduced[count] =
        (full[count] - p * reduced[count - 1]) / denominator;
    }
  } else {
    reduced[reduced_size - 1] = full[full_size - 1] / p;
    for (int count = reduced_size - 1; count >= 1; --count) {
      reduced[count - 1] =
        (full[count] - (1.0 - p) * reduced[count]) / p;
    }
  }

  return normalize_pmf(reduced);
}

inline Rcpp::NumericVector add_bernoulli_factor(
    const Rcpp::NumericVector& reduced_pmf,
    const double probability) {
  const double p = checked_probability(probability, "probability");
  Rcpp::NumericVector reduced = Rcpp::clone(reduced_pmf);
  reduced = normalize_pmf(reduced);

  const int reduced_size = reduced.size();
  Rcpp::NumericVector full(reduced_size + 1);

  for (int count = 0; count < reduced_size; ++count) {
    full[count] += (1.0 - p) * reduced[count];
    full[count + 1] += p * reduced[count];
  }

  return normalize_pmf(full);
}

inline double prior_score_from_reduced(
    const Rcpp::NumericVector& reduced_pmf,
    const double gamma) {
  if (!R_finite(gamma) || gamma <= 0.0) {
    Rcpp::stop("gamma must be positive and finite.");
  }

  double score = 0.0;
  for (R_xlen_t count = 0; count < reduced_pmf.size(); ++count) {
    score += reduced_pmf[count] *
      std::log(static_cast<double>(count) + gamma);
  }
  return score;
}

}  // namespace nsbmvi

#endif
