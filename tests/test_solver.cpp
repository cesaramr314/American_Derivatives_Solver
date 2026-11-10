// Validation tests for the amopt library.
//
// No external test framework: each test records failures through the CHECK_* macros,
// and the executable returns 1 if any check failed (this is how CTest detects it).
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "amopt/analytic.hpp"
#include "amopt/binomial.hpp"
#include "amopt/lcp.hpp"
#include "amopt/pde_solver.hpp"
#include "amopt/tridiag.hpp"

using namespace amopt;

static int g_failures = 0;

#define CHECK_TRUE(cond)                                                       \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);     \
            ++g_failures;                                                      \
        }                                                                      \
    } while (0)

#define CHECK_CLOSE(a, b, tol)                                                 \
    do {                                                                       \
        const double va_ = (a), vb_ = (b);                                     \
        if (std::abs(va_ - vb_) > (tol)) {                                     \
            std::printf("  FAIL %s:%d: %s = %.8f vs %s = %.8f (tol %.1e)\n",  \
                        __FILE__, __LINE__, #a, va_, #b, vb_, (double)(tol));  \
            ++g_failures;                                                      \
        }                                                                      \
    } while (0)

namespace {

    OptionSpec base_spec(Exercise ex) {
        OptionSpec o;
        o.S0 = 100.0;
        o.r = 0.05;
        o.q = 0.02;
        o.sigma = 0.25;
        o.T = 1.0;
        o.exercise = ex;
        return o;
    }

    GridConfig fine_grid() {
        GridConfig c;
        c.n_space = 800;
        c.n_time = 400;
        return c;
    }

    void test_thomas() {
        const std::size_t n = 50;
        Tridiag A(n);
        std::vector<double> x_true(n), d(n), x, work;
        for (std::size_t i = 0; i < n; ++i) {
            A.lower[i] = (i > 0) ? -1.0 : 0.0;
            A.upper[i] = (i + 1 < n) ? -1.5 : 0.0;
            A.diag[i] = 4.0 + 0.1 * static_cast<double>(i);
            x_true[i] = std::sin(0.3 * static_cast<double>(i));
        }
        for (std::size_t i = 0; i < n; ++i) {
            d[i] = A.diag[i] * x_true[i];
            if (i > 0) d[i] += A.lower[i] * x_true[i - 1];
            if (i + 1 < n) d[i] += A.upper[i] * x_true[i + 1];
        }
        thomas_solve(A, d, x, work);
        for (std::size_t i = 0; i < n; ++i) CHECK_CLOSE(x[i], x_true[i], 1e-12);
    }

    void test_european_vs_black_scholes() {
        const OptionSpec o = base_spec(Exercise::European);
        const GridConfig cfg = fine_grid();
        const double K = 100.0;
        CHECK_CLOSE(price_pde(o, CallPayoff(K), cfg).price, bs_call(o, K), 2e-3);
        CHECK_CLOSE(price_pde(o, PutPayoff(K), cfg).price, bs_put(o, K), 2e-3);
        CHECK_CLOSE(price_pde(o, StraddlePayoff(K), cfg).price, bs_straddle(o, K), 4e-3);
        CHECK_CLOSE(price_pde(o, CashOrNothingCall(K, 1.0), cfg).price, bs_cash_call(o, K, 1.0), 1e-4);
        CHECK_CLOSE(price_pde(o, CashOrNothingPut(K, 1.0), cfg).price, bs_cash_put(o, K, 1.0), 1e-4);
        CHECK_CLOSE(price_pde(o, AssetOrNothingCall(K), cfg).price, bs_asset_call(o, K), 5e-3);
    }

    void test_european_greeks() {
        // European call: Delta = e^{-qT} N(d1), Gamma = e^{-qT} phi(d1) / (S0 sigma sqrt(T)).
        const OptionSpec o = base_spec(Exercise::European);
        const double K = 100.0;
        const double sd = o.sigma * std::sqrt(o.T);
        const double d1 = (std::log(o.S0 / K) + (o.r - o.q + 0.5 * o.sigma * o.sigma) * o.T) / sd;
        const double delta = std::exp(-o.q * o.T) * norm_cdf(d1);
        const double gamma = std::exp(-o.q * o.T) * std::exp(-0.5 * d1 * d1) /
            (std::sqrt(2.0 * std::acos(-1.0)) * o.S0 * sd);
        const PricingResult r = price_pde(o, CallPayoff(K), fine_grid());
        CHECK_CLOSE(r.delta, delta, 1e-4);
        CHECK_CLOSE(r.gamma, gamma, 1e-4);
    }

    void test_american_call_no_dividends_equals_european() {
        OptionSpec o = base_spec(Exercise::American);
        o.q = 0.0;
        const double K = 100.0;
        PsorSolver psor;
        PenaltySolver pen;
        CHECK_CLOSE(price_pde(o, CallPayoff(K), fine_grid(), &psor).price, bs_call(o, K), 2e-3);
        CHECK_CLOSE(price_pde(o, CallPayoff(K), fine_grid(), &pen).price, bs_call(o, K), 2e-3);
    }

    void test_psor_and_penalty_agree() {
        const OptionSpec o = base_spec(Exercise::American);
        PsorSolver psor(1.5, 1e-10);
        PenaltySolver pen(1e-8);
        const PutPayoff put(100.0);
        const CashOrNothingPut dput(100.0, 1.0);
        CHECK_CLOSE(price_pde(o, put, fine_grid(), &psor).price,
            price_pde(o, put, fine_grid(), &pen).price, 1e-6);
        CHECK_CLOSE(price_pde(o, dput, fine_grid(), &psor).price,
            price_pde(o, dput, fine_grid(), &pen).price, 1e-6);
    }

    void test_american_vs_binomial() {
        const OptionSpec o = base_spec(Exercise::American);
        PenaltySolver pen;
        const double K = 100.0;
        const PutPayoff put(K);
        const CallPayoff call(K);
        const StraddlePayoff straddle(K);
        CHECK_CLOSE(price_pde(o, put, fine_grid(), &pen).price, binomial_crr(o, put, 4000), 3e-3);
        CHECK_CLOSE(price_pde(o, call, fine_grid(), &pen).price, binomial_crr(o, call, 4000), 3e-3);
        CHECK_CLOSE(price_pde(o, straddle, fine_grid(), &pen).price,
            binomial_crr(o, straddle, 4000), 5e-3);
    }

    void test_american_digital_vs_one_touch() {
        // The American cash-or-nothing option is a one-touch, which has a closed form.
        const OptionSpec o = base_spec(Exercise::American);
        PenaltySolver pen;
        PsorSolver psor;
        CHECK_CLOSE(price_pde(o, CashOrNothingCall(110.0, 1.0), fine_grid(), &pen).price,
            one_touch(o, 110.0, 1.0, true), 1e-4);
        CHECK_CLOSE(price_pde(o, CashOrNothingPut(90.0, 1.0), fine_grid(), &psor).price,
            one_touch(o, 90.0, 1.0, false), 1e-4);
        // American asset-or-nothing call: when S hits K the holder collects S = K.
        CHECK_CLOSE(price_pde(o, AssetOrNothingCall(110.0), fine_grid(), &pen).price,
            110.0 * one_touch(o, 110.0, 1.0, true), 1e-2);
    }

    void test_american_bounds() {
        // V_american >= V_european, and V >= payoff on the whole grid.
        const OptionSpec am = base_spec(Exercise::American);
        const OptionSpec eu = base_spec(Exercise::European);
        PenaltySolver pen;
        const PutPayoff put(100.0);
        const PricingResult ra = price_pde(am, put, fine_grid(), &pen);
        const PricingResult re = price_pde(eu, put, fine_grid());
        CHECK_TRUE(ra.price >= re.price);
        double worst = 0.0;
        for (std::size_t j = 0; j < ra.S.size(); ++j) worst = std::min(worst, ra.V[j] - put(ra.S[j]));
        CHECK_TRUE(worst > -1e-6);
    }

    void test_second_order_convergence() {
        // Errors with N, 2N, 4N: the ratio e(N)/e(2N) should approach 4.
        const OptionSpec o = base_spec(Exercise::European);
        const double exact = bs_put(o, 100.0);
        const PutPayoff put(100.0);
        std::vector<double> err;
        for (int N : {200, 400, 800}) {
            GridConfig c;
            c.n_space = N;
            c.n_time = N / 2;
            err.push_back(std::abs(price_pde(o, put, c).price - exact));
        }
        const double ratio1 = err[0] / err[1], ratio2 = err[1] / err[2];
        std::printf("  error ratios: %.2f, %.2f\n", ratio1, ratio2);
        CHECK_TRUE(ratio1 > 3.0 && ratio1 < 5.0);
        CHECK_TRUE(ratio2 > 3.0 && ratio2 < 5.0);
    }

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        { "thomas", test_thomas },
        { "european_vs_black_scholes", test_european_vs_black_scholes },
        { "european_greeks", test_european_greeks },
        { "american_call_no_dividends_equals_european", test_american_call_no_dividends_equals_european },
        { "psor_and_penalty_agree", test_psor_and_penalty_agree },
        { "american_vs_binomial", test_american_vs_binomial },
        { "american_digital_vs_one_touch", test_american_digital_vs_one_touch },
        { "american_bounds", test_american_bounds },
        { "second_order_convergence", test_second_order_convergence },
    };
    for (const auto& [name, fn] : tests) {
        const int before = g_failures;
        std::printf("[ RUN  ] %s\n", name.c_str());
        fn();
        std::printf("[ %s ] %s\n", g_failures == before ? " OK " : "FAIL", name.c_str());
    }
    std::printf("\n%s (%d failures)\n", g_failures == 0 ? "ALL PASSED" : "SOME TESTS FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
