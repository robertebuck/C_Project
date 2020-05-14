// cpp file for class of SDE
#include "SDE.hpp"

// constructors & destructors
SDE::SDE() : data(new OptionData(0,0,0,0,0,1)) {};  // deafult constructor
SDE::SDE(const OptionData& optionData) : data(new OptionData(optionData)) {}; 
SDE::SDE(const SDE& S) : data(new OptionData(*S.data)) {}; // copy constructor
SDE::~SDE() {}; // destructor

// Member functions

// Drift euler term
double SDE::drift_euler(double t, double S)
{ 
	return (data->r - data->D) * S; // r - D
}

// Diffusion euler term
double SDE::diffusion_euler(double t, double S)
{

	return data->sig * S;
}

// Deterministic component of exact solution
double SDE::deterministic_exact()
{
	return (data->r - data->D - 0.5*(data->sig)*(data->sig)); // r - D - 0.5*sig^2
}