// hpp file for class of SDE.
// Built upon from the SDE class given in TestMCOMP.cpp.
// Added deterministic_exact() fucntion for use in exact simulations

#ifndef SDE_HPP
#define SDE_HPP

#include "OptionData.hpp"

class SDE
{
private:

	std::shared_ptr<OptionData> data;	// The data for the option

public:

	// constructors & destructors
	SDE(); // deafult constructor
	SDE(const OptionData& optionData); 
	SDE(const SDE& S); // copy constructor
	~SDE(); // destructor
	
	// member functions
	double drift_euler(double t, double S);
	double diffusion_euler(double t, double S);
	double deterministic_exact();

};

#endif