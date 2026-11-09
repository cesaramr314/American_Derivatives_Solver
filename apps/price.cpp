// Tabla de precios para varios payoffs: europeo (BS vs EDP) y americano
// (PSOR vs penalty vs binomial).
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "amopt/analytic.hpp"
#include "amopt/binomial.hpp"
#include "amopt/lcp.hpp"
#include "amopt/pde_solver.hpp"

using namespace amopt;

int main() {
    // Binarias con strike lejos de S0: en S0 = K la americana binaria es trivial
    // (S toca K de inmediato y el valor es casi `cash`).
    const double K = 100.0, K_up = 110.0, K_dn = 90.0, cash = 1.0;
    OptionSpec opt;
    opt.S0 = 100.0;
    opt.r = 0.05;
    opt.q = 0.02;
    opt.sigma = 0.25;
    opt.T = 1.0;

    GridConfig cfg;
    cfg.n_space = 800;
    cfg.n_time = 400;

    struct Case {
        std::unique_ptr<Payoff> payoff;
        double bs;
    };
    std::vector<Case> cases;
    cases.push_back({ std::make_unique<CallPayoff>(K), bs_call(opt, K) });
    cases.push_back({ std::make_unique<PutPayoff>(K), bs_put(opt, K) });
    cases.push_back({ std::make_unique<CashOrNothingCall>(K_up, cash), bs_cash_call(opt, K_up, cash) });
    cases.push_back({ std::make_unique<CashOrNothingPut>(K_dn, cash), bs_cash_put(opt, K_dn, cash) });
    cases.push_back({ std::make_unique<AssetOrNothingCall>(K_up), bs_asset_call(opt, K_up) });
    cases.push_back({ std::make_unique<StraddlePayoff>(K), bs_straddle(opt, K) });

    PsorSolver psor(1.5, 1e-8);
    PenaltySolver penalty(1e-8);

    std::printf("S0=%.1f r=%.3f q=%.3f sigma=%.2f T=%.2f  grid %d x %d\n\n", opt.S0, opt.r,
        opt.q, opt.sigma, opt.T, cfg.n_space, cfg.n_time);
    std::printf("%-22s %5s | %9s %9s | %9s %9s %9s\n", "payoff", "K", "Eur BS", "Eur PDE",
        "Am PSOR", "Am pen.", "Am CRR");
    std::printf("%s\n", std::string(82, '-').c_str());

    for (const auto& c : cases) {
        OptionSpec eur = opt;
        eur.exercise = Exercise::European;
        OptionSpec am = opt;
        am.exercise = Exercise::American;

        const double v_eur = price_pde(eur, *c.payoff, cfg).price;
        const double v_psor = price_pde(am, *c.payoff, cfg, &psor).price;
        const double v_pen = price_pde(am, *c.payoff, cfg, &penalty).price;
        const double v_crr = binomial_crr(am, *c.payoff, 2000);

        std::printf("%-22s %5.0f | %9.5f %9.5f | %9.5f %9.5f %9.5f\n",
            c.payoff->name().c_str(), c.payoff->strike(), c.bs, v_eur, v_psor, v_pen,
            v_crr);
    }
    return 0;
}
