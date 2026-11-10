#pragma once

#include <string>
#include <vector>

#include "amopt/tridiag.hpp"

namespace amopt {
	// Solves the LCP Au - b >= 0, u - g >= 0, (Au - b)^T(u - g) = 0 in a given time step.
	// u enters with an initial guess (previous time step) and leaves with the solution.
	class LcpSolver {
	public:
		virtual ~LcpSolver() = default;
		virtual int solve(const Tridiag& A, const std::vector<double>& b,
			const std::vector<double>& g, std::vector<double>& u) = 0;
		virtual std::string name() const = 0;
	};

	// Projected SOR [Cryer, 1971]. If adaptive, omega is tuned between calls
	// (time steps) by hill climbing on the number of sweeps.
	class PsorSolver final : public LcpSolver {
	public:
		explicit PsorSolver(double omega = 1.5, double tol = 1e-8, int max_iter = 100000,
			bool adaptive = false)
			: omega_(omega), tol_(tol), max_iter_(max_iter), adaptive_(adaptive) { }
		int solve(const Tridiag& A, const std::vector<double>& b,
			const std::vector<double>& g, std::vector<double>& u) override;
		std::string name() const override { return adaptive_ ? "psor-adaptive" : "psor"; }

	private:
		double omega_, tol_;
		int max_iter_;
		bool adaptive_;
		int last_iters_ = 0;   // sweeps used in the previous call
		double step_ = 0.02;   // current change of omega (its sign is the direction)
	};

	// Penalty method [Forsyth et.al.,  2002]; penalty factor should be ~ 1/tol.
	// Keep tol >= 1e-8: larger penalties can make the active set cycle due to round-off.
	class PenaltySolver final : public LcpSolver {
	public:
		explicit PenaltySolver(double tol = 1e-8, int max_iter = 100)
			: tol_(tol), max_iter_(max_iter) {
		}
		int solve(const Tridiag& A, const std::vector<double>& b,
			const std::vector<double>& g, std::vector<double>& u) override;
		std::string name() const override { return "penalty"; }

	private:
		double tol_;
		int max_iter_;
		// Reused buffers between calls (no allocation inside the time loop)
		Tridiag M_;
		std::vector<double> rhs_, unew_, work_;
		std::vector<char> active_; // penalized nodes
	};

} // namespace amopt
