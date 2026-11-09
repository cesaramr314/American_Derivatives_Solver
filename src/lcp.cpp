#include "amopt/lcp.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace amopt {
	// Projected SOR
	int PsorSolver::solve(const Tridiag& A, const std::vector<double>& b,
		const std::vector<double>& g, std::vector<double>& u) {

		const std::size_t n = b.size();
		for (int k = 1; k <= max_iter_; k++) {
			double err = 0.0;
			for (std::size_t i = 0; i < n; i++) {
				double res = b[i] - A.diag[i] * u[i];
				if (i > 0) res -= A.lower[i] * u[i - 1];
				if (i + 1 < n) res -= A.upper[i] * u[i + 1];

				const double unew = std::max(g[i], u[i] + omega_ * res / A.diag[i]);
				err = std::max(err, std::abs(unew - u[i]) / std::max(1.0, std::abs(unew)));
				u[i] = unew;
			}
			if (err < tol_) return k;
		}
		throw std::runtime_error("PSOR did not converge. Increase max_iter or check omega");
	}


	// Penalty
	int PenaltySolver::solve(const Tridiag& A, const std::vector<double>& b,
		const std::vector<double>& g, std::vector<double>& u) {

		const std::size_t n = b.size();
		const double big = 1.0 / tol_;
		M_ = A;
		rhs_.resize(n);

		active_.assign(n, 0);

		// Penalize the nodes where the current iterate violates u >= g
		for (int k = 1; k <= max_iter_; k++) {
			bool same_set = (k > 1);
			for (std::size_t i = 0; i < n; i++) {
				const char act = (u[i] < g[i]) ? 1 : 0;
				if (act != active_[i]) same_set = false;
				active_[i] = act;
				const double p = act ? big : 0.0;
				M_.diag[i] = A.diag[i] + p;
				rhs_[i] = b[i] + p * g[i];
			}
			// If the penalized set did not change, we have convergence
			if (same_set) return k - 1;

			thomas_solve(M_, rhs_, unew_, work_);

			double err = 0;

			for (std::size_t i = 0; i < n; i++) {
				err = std::max(err, std::abs((unew_[i] - u[i]) / std::max(1.0, std::abs(unew_[i]))));
			}
			u.swap(unew_);
			if (err < tol_) return k;
		}
		throw std::runtime_error("Penalty did not converge. Increase max_iter or check M-matrix condition");
	}

} // namespace amopt