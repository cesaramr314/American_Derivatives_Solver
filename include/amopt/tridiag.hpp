#pragma once

#include <cstddef>
#include <vector>

namespace amopt {
	// Stores the main diagonals of a tridiagonal matrix in a compact form.
	// Row i: lower[i]*x[i] + diag[i]*x[i] + upper[i]*x[i+1]
	// lower[0] and upper[n-1] are not used.

	struct Tridiag {
		std::vector <double> lower, diag, upper;

		explicit Tridiag(std::size_t n=0) : lower(n,0.0), diag(n,0.0), upper(n,0.0) {}
		std::size_t size() const { return diag.size(); }
	};

	// Solves a tridiagonal system Ax = d using the Thomas algorithm in O(n).
	// cp is a buffer variable.
	// Requires that A does not need pivoting (i.e., no zero diagonal entries).

	void thomas_solve(const Tridiag& A, const std::vector<double>& d,
		std::vector<double>& x, std::vector<double>& cp);

}