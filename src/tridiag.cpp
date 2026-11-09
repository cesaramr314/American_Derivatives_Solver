#include "amopt/tridiag.hpp"

namespace amopt {

    void thomas_solve(const Tridiag& A, const std::vector<double>& d,
        std::vector<double>& x, std::vector<double>& cp) {
        const std::size_t n = A.size();
        x.resize(n);
        cp.resize(n);

		// Forward elimination.
        double beta = A.diag[0];
        x[0] = d[0] / beta;
        for (std::size_t i = 1; i < n; i++) {
            cp[i - 1] = A.upper[i - 1] / beta;
            beta = A.diag[i] - A.lower[i] * cp[i - 1];
            x[i] = (d[i] - A.lower[i] * x[i - 1]) / beta;
        }
		// Back substitution.
        for (std::size_t i = n - 1; i > 0; i--) {
            x[i - 1] -= cp[i - 1] * x[i];
        }
    }

}  
