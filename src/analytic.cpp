#include "amopt/analytic.hpp"

#include <cmath>

namespace amopt {

    namespace {
        struct D12 {
            double d1, d2;
        };
        D12 d12(const OptionSpec& o, double K) {
            const double sd = o.sigma * std::sqrt(o.T);
            const double d1 = (std::log(o.S0 / K) + (o.r - o.q + 0.5 * o.sigma * o.sigma) * o.T) / sd;
            return { d1, d1 - sd };
        }
    }  // namespace

    double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

    double bs_call(const OptionSpec& o, double K) {
        const auto [d1, d2] = d12(o, K);
        return o.S0 * std::exp(-o.q * o.T) * norm_cdf(d1) - K * std::exp(-o.r * o.T) * norm_cdf(d2);
    }

    double bs_put(const OptionSpec& o, double K) {
        const auto [d1, d2] = d12(o, K);
        return K * std::exp(-o.r * o.T) * norm_cdf(-d2) - o.S0 * std::exp(-o.q * o.T) * norm_cdf(-d1);
    }

    double bs_cash_call(const OptionSpec& o, double K, double cash) {
        return cash * std::exp(-o.r * o.T) * norm_cdf(d12(o, K).d2);
    }

    double bs_cash_put(const OptionSpec& o, double K, double cash) {
        return cash * std::exp(-o.r * o.T) * norm_cdf(-d12(o, K).d2);
    }

    double bs_asset_call(const OptionSpec& o, double K) {
        return o.S0 * std::exp(-o.q * o.T) * norm_cdf(d12(o, K).d1);
    }

    double bs_straddle(const OptionSpec& o, double K) { return bs_call(o, K) + bs_put(o, K); }

    double one_touch(const OptionSpec& o, double K, double cash, bool is_call) {
        const double s2 = o.sigma * o.sigma;
        const double sd = o.sigma * std::sqrt(o.T);
        const double mu = (o.r - o.q - 0.5 * s2) / s2;
        const double lam = std::sqrt(mu * mu + 2.0 * o.r / s2);
        const double h = K / o.S0;
        const double z = std::log(h) / sd + lam * sd;
        const double eta = is_call ? -1.0 : 1.0;  // barrier up (call) or down (put)
        return cash * (std::pow(h, mu + lam) * norm_cdf(eta * z) +
            std::pow(h, mu - lam) * norm_cdf(eta * (z - 2.0 * lam * sd)));
    }

}  // namespace amopt