#pragma once

#include <algorithm>
#include <cmath>
#include <string>

namespace amopt {

    // Common interface for option payoffs.
    //  - operator()(S): intrinsic value g(S).
    //  - far_field(S, tau, r, q): asymptotic value of the option as S -> 0 or S -> infinity,
    //    used as Dirichlet boundary conditions in the PDE solver.
    class Payoff {
    public:
        virtual ~Payoff() = default;
        virtual double operator()(double S) const = 0;
        virtual double far_field(double S, double tau, double r, double q) const = 0;
        virtual double strike() const = 0;
        virtual bool is_discontinuous() const { return false; }
        virtual std::string name() const = 0;
    };

    class CallPayoff final : public Payoff {
    public:
        explicit CallPayoff(double K) : K_(K) {}
        double operator()(double S) const override { return std::max(S - K_, 0.0); }
        double far_field(double S, double tau, double r, double q) const override {
            return std::max(S * std::exp(-q * tau) - K_ * std::exp(-r * tau), 0.0);
        }
        double strike() const override { return K_; }
        std::string name() const override { return "call"; }

    private:
        double K_;
    };

    class PutPayoff final : public Payoff {
    public:
        explicit PutPayoff(double K) : K_(K) {}
        double operator()(double S) const override { return std::max(K_ - S, 0.0); }
        double far_field(double S, double tau, double r, double q) const override {
            return std::max(K_ * std::exp(-r * tau) - S * std::exp(-q * tau), 0.0);
        }
        double strike() const override { return K_; }
        std::string name() const override { return "put"; }

    private:
        double K_;
    };

    // Pays `cash` if S > K.
    class CashOrNothingCall final : public Payoff {
    public:
        CashOrNothingCall(double K, double cash) : K_(K), cash_(cash) {}
        double operator()(double S) const override { return S > K_ ? cash_ : 0.0; }
        double far_field(double S, double tau, double r, double /*q*/) const override {
            return S > K_ ? cash_ * std::exp(-r * tau) : 0.0;
        }
        double strike() const override { return K_; }
        bool is_discontinuous() const override { return true; }
        std::string name() const override { return "cash-or-nothing call"; }

    private:
        double K_, cash_;
    };

    // Pays `cash` if S < K.
    class CashOrNothingPut final : public Payoff {
    public:
        CashOrNothingPut(double K, double cash) : K_(K), cash_(cash) {}
        double operator()(double S) const override { return S < K_ ? cash_ : 0.0; }
        double far_field(double S, double tau, double r, double /*q*/) const override {
            return S < K_ ? cash_ * std::exp(-r * tau) : 0.0;
        }
        double strike() const override { return K_; }
        bool is_discontinuous() const override { return true; }
        std::string name() const override { return "cash-or-nothing put"; }

    private:
        double K_, cash_;
    };

    // Pays S if S > K.
    class AssetOrNothingCall final : public Payoff {
    public:
        explicit AssetOrNothingCall(double K) : K_(K) {}
        double operator()(double S) const override { return S > K_ ? S : 0.0; }
        double far_field(double S, double tau, double /*r*/, double q) const override {
            return S > K_ ? S * std::exp(-q * tau) : 0.0;
        }
        double strike() const override { return K_; }
        bool is_discontinuous() const override { return true; }
        std::string name() const override { return "asset-or-nothing call"; }

    private:
        double K_;
    };

    // |S - K| = call + put with the same strike.
    class StraddlePayoff final : public Payoff {
    public:
        explicit StraddlePayoff(double K) : K_(K) {}
        double operator()(double S) const override { return std::abs(S - K_); }
        double far_field(double S, double tau, double r, double q) const override {
            return std::abs(S * std::exp(-q * tau) - K_ * std::exp(-r * tau));
        }
        double strike() const override { return K_; }
        std::string name() const override { return "straddle"; }

    private:
        double K_;
    };

}  // namespace amopt