// Compares PSOR (fixed and adaptive omega) and the penalty method on an American put:
// error, iterations and time as the grid is refined (dt proportional to dx).
// Writes a CSV to the file given as argument (or to the screen if none is given).
//
//   ./amopt_benchmark results/benchmark.csv
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "amopt/lcp.hpp"
#include "amopt/pde_solver.hpp"

using namespace amopt;

int main(int argc, char* argv[]) {
    // Output: the file given as argument, or the screen.
    std::ofstream file;
    std::ostream* out = &std::cout;
    if (argc > 1) {
        file.open(argv[1]);
        if (!file) {
            std::cerr << "could not open " << argv[1] << " (does the folder exist?)\n";
            return 1;
        }
        out = &file;
    }
    char line[256];

    OptionSpec opt;
    opt.S0 = 100.0;
    opt.r = 0.05;
    opt.q = 0.0;
    opt.sigma = 0.3;
    opt.T = 1.0;
    opt.exercise = Exercise::American;
    const PutPayoff put(100.0);

    // Reference: penalty method on a very fine grid.
    PenaltySolver ref_solver(1e-8);
    GridConfig ref_cfg;
    ref_cfg.n_space = 16000;
    ref_cfg.n_time = 8000;
    const double ref = price_pde(opt, put, ref_cfg, &ref_solver).price;
    std::fprintf(stderr, "reference (penalty, %d x %d) = %.10f\n", ref_cfg.n_space,
        ref_cfg.n_time, ref);

    *out << "method,n_space,n_time,price,abs_error,iters_per_step,seconds\n";
    for (int N : {100, 200, 400, 800, 1600, 3200}) {
        GridConfig cfg;
        cfg.n_space = N;
        cfg.n_time = N / 2;

        for (int method = 0; method < 3; ++method) {
            // Best of 3 runs, with fresh solvers in each run so that the adaptive
            // omega does not carry over what it learned in the previous run.
            PricingResult best;
            best.seconds = 1e300;
            std::string name;
            for (int rep = 0; rep < 3; ++rep) {
                PsorSolver psor(1.5, 1e-8);
                PsorSolver psor_adaptive(1.5, 1e-8, 100000, true);
                PenaltySolver penalty(1e-8);
                LcpSolver* solvers[] = { &psor, &psor_adaptive, &penalty };
                LcpSolver* solver = solvers[method];
                name = solver->name();
                PricingResult r = price_pde(opt, put, cfg, solver);
                if (r.seconds < best.seconds) best = r;
            }
            const int steps = cfg.n_time + cfg.rannacher_steps;  // each Rannacher step is 2 half steps
            std::snprintf(line, sizeof(line), "%s,%d,%d,%.10f,%.3e,%.2f,%.6f\n",
                name.c_str(), N, cfg.n_time, best.price,
                std::abs(best.price - ref),
                static_cast<double>(best.iterations) / steps, best.seconds);
            *out << line;
        }
    }
    return 0;
}