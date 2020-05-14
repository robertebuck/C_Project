// Builds upon OptionData.hpp file given in class.
// Addition of two new classes and two new payoffs to the function myPayOffFunction


#ifndef OptionData_HPP
#define OptionData_HPP

#include <iostream>
#include <algorithm> // for max()
#include <boost/parameter.hpp>
#include <vector>
#include <numeric>

namespace OptionParams
{
	BOOST_PARAMETER_KEYWORD(Tag, strike)
	BOOST_PARAMETER_KEYWORD(Tag, expiration)
	BOOST_PARAMETER_KEYWORD(Tag, interestRate)
	BOOST_PARAMETER_KEYWORD(Tag, volatility)
	BOOST_PARAMETER_KEYWORD(Tag, dividend)
	BOOST_PARAMETER_KEYWORD(Tag, optionType)
}


// Encapsulate all data in one place
struct OptionData
{

	double K; // strike
	double T; // time to expiration
	double r; // risk free rate
	double sig; // annual volatility
	int type; // option type: 1 == european call, 2 == european put, 3 == arithmetic asian call, 4 == arithmetic asian put, deafult european call

	// Extra data 
	double D; // dividend

	// constructor
	explicit constexpr OptionData(double strike, double expiration, double interestRate,
									double volatility, double dividend, int PC)
		: K(strike), T(expiration), r(interestRate), sig(volatility), D(dividend), type(PC)
	{}

	//copy constructor
	explicit constexpr OptionData(const OptionData& copy)
		: K(copy.K), T(copy.T), r(copy.r), sig(copy.sig), D(copy.D), type(copy.type)
	{}
	
	template <typename ArgPack> OptionData(const ArgPack& args)
	{
		K = args[OptionParams::strike];
		T = args[OptionParams::expiration];
		r = args[OptionParams::interestRate];
		sig = args[OptionParams::volatility];
		D = args[OptionParams::dividend];
		type = args[OptionParams::optionType];
	}

	// Payoff function
	double myPayOffFunction(std::vector<double>& S)
	{

		// European Put
		if (type == 2)
		{
			return std::max(K - S.back(), 0.0);
		}
		
		// Arithmetic asian call
		else if (type == 3)
		{
			double average = accumulate(S.begin(), S.end(), 0.0) / S.size();
			
			return std::max(average - K, 0.0);
		}

		// Arithmetic asian put
		else if (type == 4)
		{
			double average = accumulate(S.begin(), S.end(), 0.0) / S.size();

			return std::max(K - average, 0.0);
		}

		// European Call
		else
		{ 
			return std::max(S.back() - K, 0.0);
		}
	}
};


#endif