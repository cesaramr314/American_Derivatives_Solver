#include "amopt/pde_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "amopt/tridiag.hpp"

namespace amopt {
	PricingResult price_pde(const OptionSpec& opt, const Payoff& payoff,
		const GridConfig& cfg, LcpSolver* lcp) {

		if (cfg.n_space < 4 || cfg.n_space % 2 != 0) {
			throw std::invalid_argument("n_space must be even and >= 4");
		}
		if (cfg.n_time < 1) {
			throw std::invalid_argument("n_time must be positive");
		}
		const bool american = (opt.exercise == Exercise::American);
		if (american && lcp == nullptr) {
			throw std::invalid_argument("lcp solver must be provided for American options");
		}

		const auto t_start = std::chrono::steady_clock::now();


		// ---- Mesh in x = ln S centered at K (K is in node N/2)
		const int N = cfg.n_space;
		const std::size_t n = static_cast<std::size_t>(N) + 1;
		const double K = payoff.strike();
		const double x0 = std::log(opt.S0);
		const double xK = std::log(K);
		const double sd = opt.sigma * std::sqrt(opt.T);
		const double L = std::max(cfg.width * sd, std::abs(x0 - xK) + 2.0 * sd);
		const double dx = 2.0 * L / N;
		const double xmin = xK - L;

		// g: obstacle restricts u >= g. u: initial condition at tau = 0
		std::vector<double> S(n), g(n);
		for (std::size_t j = 0; j < n; j++) {
			S[j] = std::exp(xmin + static_cast<double>(j) * dx);
			g[j] = payoff(S[j]);
		}
		std::vector<double> u = g;

		// Discontinuous payoff at K (node N/2): the initial condition takes the average of
		// the one-sided limits [Pooley et.al.,  2003]; the obstacle u >= g
		// [Wilmott et.al., 1993] takes their maximum, since at S = K the holder
		// can exercise and collect the higher side.
		if (payoff.is_discontinuous()) {
			const std::size_t k = static_cast<std::size_t>(N / 2);
			const double left = payoff(K * (1.0 - 1e-12));
			const double right = payoff(K * (1.0 + 1e-12));
			u[k] = 0.5 * (left + right);
			g[k] = std::max(left, right);
		}



		// ---- Difference scheme with central differences; L u = a u_xx + b u_x - r u.
		const double a = 0.5 * opt.sigma * opt.sigma;
		const double b = opt.r - opt.q - a;
		const double lo = a / (dx * dx) - b / (2.0 * dx);
		const double di = -2.0 * a / (dx * dx) - opt.r;
		const double up = a / (dx * dx) + b / (2.0 * dx);

		// Left hand side matrix of the scheme  I - theta*h*L  (boundary rows are identities).
		auto make_lhs = [&](double theta, double h) {
			Tridiag A(n);
			A.diag[0] = 1.0;
			A.diag[n - 1] = 1.0;
			for (std::size_t i = 1; i < n - 1; i++) {
				A.lower[i] = -theta * h * lo;
				A.diag[i] = 1.0 - theta * h * di;
				A.upper[i] = -theta * h * up;
			}
			return A;
			};

		const double dt = opt.T / cfg.n_time;
		const Tridiag A_cn = make_lhs(0.5, dt); // Crank-Nicolson
		const Tridiag A_half = make_lhs(1.0, 0.5 * dt); // half step of fully implicit



		// ---- Time stepping on tau = T - t
		std::vector<double> rhs(n), work(n);
		double tau = 0.0;
		long iterations = 0;


		auto step = [&](const Tridiag& A, double theta, double h) {
			tau += h;
			const double e = (1.0 - theta) * h;
			for (std::size_t i = 1; i < n - 1; i++) {
				rhs[i] = u[i] + e * (lo * u[i - 1] + di * u[i] + up * u[i + 1]);
			}
			rhs[0] = payoff.far_field(S[0], tau, opt.r, opt.q);
			rhs[n - 1] = payoff.far_field(S[n - 1], tau, opt.r, opt.q);
			if (american) {
				rhs[0] = std::max(rhs[0], g[0]);
				rhs[n - 1] = std::max(rhs[n - 1], g[n - 1]);
				iterations += lcp->solve(A, rhs, g, u);
			}
			else {
				thomas_solve(A, rhs, u, work);
				iterations += 1;
			}
			};

		for (int m = 0; m < cfg.n_time; m++) {
			if (m < cfg.rannacher_steps) {
				step(A_half, 1.0, 0.5 * dt);
				step(A_half, 1.0, 0.5 * dt);
			}
			else {
				step(A_cn, 0.5, dt);
			}
		}


		// ---- Quadratic interpolation in x0 and greeks ----
		long ic = std::lround((x0 - xmin) / dx);
		ic = std::clamp(ic, 1L, static_cast<long>(N) - 1);
		const std::size_t i = static_cast<std::size_t>(ic);
		const double t = (x0 - (xmin + static_cast<double>(i) * dx)) / dx;
		const double d1 = 0.5 * (u[i + 1] - u[i - 1]);
		const double d2 = u[i + 1] - 2.0 * u[i] + u[i - 1];
		const double ux = (d1 + t * d2) / dx;
		const double uxx = d2 / (dx * dx);

		// V(S) = u(ln S), so we need to use chain rule for the greeks
		PricingResult res;
		res.price = u[i] + t * d1 + 0.5 * t * t * d2;
		res.delta = ux / opt.S0;
		res.gamma = (uxx - ux) / (opt.S0 * opt.S0);
		res.iterations = iterations;
		res.m_matrix = (lo >= 0.0 && up >= 0.0);
		res.S = std::move(S);
		res.V = std::move(u);
		res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
		return res;


	}

}	// namespace amopt