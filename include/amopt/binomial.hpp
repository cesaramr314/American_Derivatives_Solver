#pragma once

#include "amopt/option.hpp"
#include "amopt/payoff.hpp"

namespace amopt {

	// Árbol binomial de Cox-Ross-Rubinstein. Promedia n y n+1 pasos para
	// amortiguar la oscilación par/impar. Sirve como referencia independiente.
	double binomial_crr(const OptionSpec& o, const Payoff& payoff, int steps);

}  // namespace amopt