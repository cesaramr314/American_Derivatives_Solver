#pragma once

#include "amopt/option.hpp"
#include "amopt/payoff.hpp"

namespace amopt {

	// Cox-Ross-Rubinstein binomial tree. Averages n and n+1 steps to damp
	// the even/odd oscillation. Used as an independent reference.
	double binomial_crr(const OptionSpec& o, const Payoff& payoff, int steps);

}  // namespace amopt
