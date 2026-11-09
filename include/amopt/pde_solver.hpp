#pragma once

#include <vector>

#include "amopt/lcp.hpp"
#include "amopt/option.hpp"
#include "amopt/payoff.hpp"

namespace amopt {

    struct GridConfig {
        int n_space = 400;        // number of space intervals in x = ln S (must be even)
        int n_time = 200;         // number of time steps
        double width = 6.0;       // half width of the grid in standard deviations sigma*sqrt(T)
        int rannacher_steps = 2;  // CN steps replaced by two implicit half steps each
    };

    struct PricingResult {
        double price = 0.0;
        double delta = 0.0;
        double gamma = 0.0;
        long iterations = 0;    // total LCP solver iterations (summed over time steps)
        double seconds = 0.0;   // wall time
        bool m_matrix = true;   // whether the discretization is an M-matrix
        std::vector<double> S;  // grid in S
        std::vector<double> V;  // values at tau = T
    };

    // Prices the option by solving the Black-Scholes PDE in x = ln S with
    // Crank-Nicolson + Rannacher start. American exercise uses lcp at every
    // step. For European exercise it is ignored and may be nullptr.
    PricingResult price_pde(const OptionSpec& opt, const Payoff& payoff,
        const GridConfig& cfg, LcpSolver* lcp = nullptr);

}  // namespace amopt