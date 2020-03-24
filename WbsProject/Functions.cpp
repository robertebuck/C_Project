// Cpp file for Functions header file

#include "Functions.hpp"

    std::default_random_engine dre;
    std::normal_distribution<double> nor(0.0, 1.0);
void generateRN(std::vector<double>& v)
{
    // loop to generate vector of standard normal random numbers
	for (std::size_t j = 0; j < v.size(); ++j)
	{
		v[j] = nor(dre);
	}
}

// Normal variates
double n(double x)
{
	double A = 1.0 / std::sqrt(2.0 * 3.1415);
	return A * std::exp(-x * x * 0.5);
}

double N(double x)
{ // The approximation to the cumulative normal distribution

	return 0.5 * (1.0 - std::erf(-x / std::sqrt(2.0)));
}

void data_output_euler(std::vector<std::vector<double>> temp, long NSIM, long NT)
{ //Creating excel file for paths from euler simulation
    std::ofstream out("paths_euler.csv");
    // looping through paths to output in excel
    for (int j = 0; j <= NT; j++)
    {
        for (int i = 0; i < NSIM; i++)
        {
            out << temp[i][j] << ',';
        }
        out << '\n';
    }
}

void data_output_exact1(std::vector<std::vector<double>> temp, long NSIM, long NT)
{ //Creating excel file for paths from exact simulation
    std::ofstream out("paths_exact1.csv");
    // looping through paths to output in excel
    for (int j = 0; j <= NT; j++)
    {
        for (int i = 0; i < NSIM; i++)
        {
            out << temp[i][j] << ',';
        }
        out << '\n';
    }
}


void data_output_exact2(std::vector<std::vector<double>> temp, long NSIM, long NT)
{ //Creating excel file for paths from exact simulation
    std::ofstream out("paths_exact2.csv");
    // looping through paths to output in excel
    for (int j = 0; j <= NT; j++)
    {
        for (int i = 0; i < NSIM; i++)
        {
            out << temp[i][j] << ',';
        }
        out << '\n';
    }
}


