#pragma once

#include "amopt/option.hpp"

namespace amopt {

	// Black-Scholes closed-form solutions for European options and related derivatives.
	double norm_cdf(double x);
	double bs_call(const OptionSpec& o, double K);
	double bs_put(const OptionSpec& o, double K);
	double bs_cash_call(const OptionSpec& o, double K, double cash);
	double bs_cash_put(const OptionSpec& o, double K, double cash);
	double bs_asset_call(const OptionSpec& o, double K);
	double bs_straddle(const OptionSpec& o, double K);

	// American Binary cash-or-nothing ("one-touch"). Closed form from [Reiner-Rubinstein, 1991].
	double one_touch(const OptionSpec& o, double K, double cash, bool is_call);

}  // namespace amopt