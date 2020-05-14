// Header file of functions used in the project

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <random>
#include <vector>
#include <fstream>

// Generates a vector of standard normally distributed random numbers
void generateRN(std::vector<double>& v);

// Normal variates
double n(double x);

// The approximation to the cumulative normal distribution
double N(double x);

// Outputs a CSV file for a 2D vector of euler simulations
void data_output_euler(std::vector<std::vector<double>> temp, long NSIM, long NT);

// Outputs a CSV file for a 2D vector of exact1 simulations
void data_output_exact1(std::vector<std::vector<double>> temp, long NSIM, long NT);

// Outputs a CSV file for a 2D vector of exact simulations
void data_output_exact2(std::vector<std::vector<double>> temp, long NSIM, long NT);

// Output a CSV file for a range of option prices for a range of spot prices using the euler simulation
void data_output_euler_range(std::vector<std::vector<double>> temp, std::vector<double> S, int f);

// Output a CSV file for a range of option prices for a range of spot prices using the Exact1 simluation
void data_output_exact1_range(std::vector<std::vector<double>> temp, std::vector<double> S, int f);

// Output a CSV file for a range of option prices for a range of spot prices using the Exact2 simluation
void data_output_exact2_range(std::vector<std::vector<double>> temp, std::vector<double> S, int f);

#endif

