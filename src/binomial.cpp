#include "amopt/binomial.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace amopt {

    namespace {
        double crr_once(const OptionSpec& o, const Payoff& payoff, int n) {
            const double dt = o.T / n;
            const double s = o.sigma * std::sqrt(dt);
            const double up = std::exp(s), dn = 1.0 / up;
            const double p = (std::exp((o.r - o.q) * dt) - dn) / (up - dn);
            const double disc = std::exp(-o.r * dt);
            const bool american = (o.exercise == Exercise::American);

            std::vector<double> v(static_cast<std::size_t>(n) + 1);
            for (int j = 0; j <= n; ++j) v[j] = payoff(o.S0 * std::exp((2 * j - n) * s));

            for (int step = n - 1; step >= 0; --step) {
                for (int j = 0; j <= step; ++j) {
                    v[j] = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
                    if (american) v[j] = std::max(v[j], payoff(o.S0 * std::exp((2 * j - step) * s)));
                }
            }
            return v[0];
        }
    }  // namespace

    double binomial_crr(const OptionSpec& o, const Payoff& payoff, int steps) {
        return 0.5 * (crr_once(o, payoff, steps) + crr_once(o, payoff, steps + 1));
    }

}  // namespace amopt