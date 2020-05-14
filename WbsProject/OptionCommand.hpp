// OptionCommand.hpp
// Taken directly from OptionCommand.hpp file given in class.
// Small modifications made moving the N and n fucntions to Functions.hpp and Functions.cpp
// and making the output from each execute command a double so that comparision  of simulated price
// for European options can be made easier.
// Note not seperated into hpp and cpp file as use of Black Scholes formula to compute European
// option prices is outside the scope of Section D of the project.

#ifndef OptionCommand_HPP
#define OptionCommand_HPP

#include <memory>
#include <cmath>
#include <iostream>
#include "Functions.hpp"
#include <fstream>
#include <iterator>

class OptionCommand
{
private:

protected:
	double K;	double T;	double r;	double b;	double sig;

public:
	OptionCommand() = default;

	explicit OptionCommand(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: K(strike), T(expiration), r(riskFree), b(costOfCarry), sig(volatility) {}
	// Want to forbid copy constructor and assignment operator
	OptionCommand(const OptionCommand& c) = default;
	OptionCommand& operator = (const OptionCommand& c) = default;

	// The abstract interface
	virtual double execute(double S) = 0;

	// Implement as a function object; simple example of Template Method Pattern
	virtual double operator () (double S)
	{
		// Call derived class' execute()
		return execute(S);
	}
};

class CallPrice final : public OptionCommand
{
public:
	explicit CallPrice(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (r - b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double blacks = (S * std::exp(-b * T) * N(d1)) - (K * std::exp(-r * T) * N(d2));

		return blacks;
	}
};


class PutPrice final : public OptionCommand
{
public:
	explicit PutPrice(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (r - b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double blacks = (K * std::exp(-r * T) * N(-d2)) - (S * std::exp(-b * T) * N(-d1));

		return blacks;
	}
};

class CallDelta final : public OptionCommand
{
public:
	explicit CallDelta(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S)  override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double call_delta = std::exp((b - r) * T) * N(d1);

		return call_delta;
	}
};

class PutDelta final : public OptionCommand
{
public:
	explicit PutDelta(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double put_delta = std::exp((b - r) * T) * (N(d1) - 1.0);

		return put_delta;
	}
};

class CallGamma final : public OptionCommand
{
public:
	explicit CallGamma(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double call_gamma = n(d1) * std::exp((b - r) * T) / (S * tmp);

		return call_gamma;
	}
};


class PutGamma final : public OptionCommand
{
public:
	explicit PutGamma(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double put_gamma = n(d1) * std::exp((b - r) * T) / (S * tmp);

		return put_gamma;
	}
};

class CallVega final : public OptionCommand
{
public:
	explicit CallVega(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double call_vega = S * std::exp((b - r) * T) * n(d1) * sqrt(T);

		return call_vega;
	}
};


class PutVega final : public OptionCommand
{
public:
	explicit PutVega(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double put_vega = S * std::exp((b - r) * T) * n(d1) * sqrt(T);

		return put_vega;
	}
};

class CallRho final : public OptionCommand
{
public:
	explicit CallRho(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double call_rho = T * K * std::exp(-r * T) * N(d2);

		return call_rho;
	}
};

class PutRho final : public OptionCommand
{
public:
	explicit PutRho(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;
		double put_rho = -T * K * std::exp(-r * T) * N(-d2);

		return put_rho;
	}
};

class CallTheta final : public OptionCommand
{
public:
	explicit CallTheta(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;

		double t1 = (S * std::exp((b - r) * T) * n(d1) * sig * 0.5) / std::sqrt(T);
		double t2 = (b - r) * (S * std::exp((b - r) * T) * N(d1));
		double t3 = r * K * std::exp(-r * T) * N(d2);
		double call_theta = -(t1 + t2 + t3);

		return call_theta;
	}
};

class PutTheta final : public OptionCommand
{
public:
	explicit PutTheta(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;
		double d2 = d1 - tmp;

		double t1 = -(S * std::exp((b - r) * T) * n(d1) * sig * 0.5) / std::sqrt(T);
		double t2 = (b - r) * (S * std::exp((b - r) * T) * N(-d1));
		double t3 = r * K * std::exp(-r * T) * N(-d2);
		double put_theta = t1 + t2 + t3;

		return put_theta;
	}
};


class CallElasticity final : public OptionCommand
{
public:
	explicit CallElasticity(double strike, double expiration, double riskFree, double costOfCarry, double volatility)
		: OptionCommand(strike, expiration, riskFree, costOfCarry, volatility) {}

	virtual double execute(double S) override
	{
		compute(S, 0.25);
	}

	virtual double compute(double S, double percentageMovement)
	{
		double tmp = sig * std::sqrt(T);
		double d1 = (log(S / K) + (b + (sig * sig) * 0.5) * T) / tmp;

		double cd = std::exp((b - r) * T) * N(d1);
		double call_elasticity = (cd * S) / percentageMovement;
		
		return call_elasticity;

	}
};


#endif
