// Header file of the Class MonteCarlo

#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <iostream>
#include <memory>
#include <cmath>
#include "OptionData.hpp"
#include "SDE.hpp"
#include "Stopwatch.cpp"
#include <omp.h>
#include "Functions.hpp"


class MonteCarlo
{
private:
	std::shared_ptr<OptionData> myOption; //Info on the option being priced
	long NSIM; // Number of simulation
	long NT; // Number of time steps
	double S_0; // inital value of the underlying

public:

	// Constructors & Destructors 
	MonteCarlo(); // deafult constructor
	MonteCarlo(const OptionData& optionData, long sim, long Time, double Spot);
	MonteCarlo(const MonteCarlo& Mc); // copy constructor
	~MonteCarlo(); // detructor

	// Functions
	std::vector<double> Euler();
	std::vector<double> Exact1();
	std::vector<double> Exact2();
	void Euler_excel();
	void Exact1_excel();
	void Exact2_excel();
};

#endif

