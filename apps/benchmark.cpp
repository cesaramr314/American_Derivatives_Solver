// Compara PSOR y penalty para un put americano: error, iteraciones y tiempo
// al refinar la malla (dt proporcional a dx). Escribe un CSV en el archivo
// que se pase como argumento (o en la pantalla si no se pasa ninguno).
//
//   ./amopt_benchmark results/benchmark.csv
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "amopt/lcp.hpp"
#include "amopt/pde_solver.hpp"

using namespace amopt;

int main(int argc, char* argv[]) {
    // Salida: el archivo del argumento, o la pantalla.
    std::ofstream file;
    std::ostream* out = &std::cout;
    if (argc > 1) {
        file.open(argv[1]);
        if (!file) {
            std::cerr << "no se pudo abrir " << argv[1] << " (existe la carpeta?)\n";
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

    // Referencia: penalty en una malla muy fina.
    PenaltySolver ref_solver(1e-8);
    GridConfig ref_cfg;
    ref_cfg.n_space = 16000;
    ref_cfg.n_time = 8000;
    const double ref = price_pde(opt, put, ref_cfg, &ref_solver).price;
    std::fprintf(stderr, "referencia (penalty, %d x %d) = %.10f\n", ref_cfg.n_space,
        ref_cfg.n_time, ref);

    *out << "method,n_space,n_time,price,abs_error,iters_per_step,seconds\n";
    for (int N : {100, 200, 400, 800, 1600, 3200}) {
        GridConfig cfg;
        cfg.n_space = N;
        cfg.n_time = N / 2;

        PsorSolver psor(1.5, 1e-8);
        PenaltySolver penalty(1e-8);
        for (LcpSolver* solver : std::vector<LcpSolver*>{ &psor, &penalty }) {
            // Mejor de 3 corridas para reducir ruido de tiempo.
            PricingResult best;
            best.seconds = 1e300;
            for (int rep = 0; rep < 3; ++rep) {
                PricingResult r = price_pde(opt, put, cfg, solver);
                if (r.seconds < best.seconds) best = r;
            }
            const int steps = cfg.n_time + cfg.rannacher_steps;  // cada paso Rannacher son 2 medios pasos
            std::snprintf(line, sizeof(line), "%s,%d,%d,%.10f,%.3e,%.2f,%.6f\n",
                solver->name().c_str(), N, cfg.n_time, best.price,
                std::abs(best.price - ref),
                static_cast<double>(best.iterations) / steps, best.seconds);
            *out << line;
        }
    }
    return 0;
}