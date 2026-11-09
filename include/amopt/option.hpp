#pragma once

namespace amopt {
	enum class Exercise { European, American};

	// Market parameters

	// Options parameters. Strike price is implemented in payoff
	struct OptionSpec {
		double S0 = 100.0; // Initial stock price
		double r = 0.05; // Risk-free interest rate
		double q = 0.0; // Dividend rate
		double sigma = 0.2; // Volatility
		double T = 1.0; // Time to maturity in years
		Exercise exercise = Exercise::European; // Option type
	};
}