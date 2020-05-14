#include "MonteCarlo.hpp"
#include "OptionCommand.hpp"
#include <iomanip>

int main()
{
	double S_start;
	double S_end;
	int grid_size;
	double k;
	double T;
	double sig;
	double d;
	int type;
	double r;
	long NSIM;
	long NT;
	int s;

	std::cout << "Imput option parameters \n";
	std::cout << "Strike: ";
	std::cin >> k;
	std::cout << "Time to expiration in years: ";
	std::cin >> T;
	std::cout << "Annual volatility, e.g. 0.1 = 10%: ";
	std::cin >> sig;
	std::cout << "Annual cost of carry/yearly dividend, e.g. 0.1 = 10%: ";
	std::cin >> d;
	std::cout << "Annual risk free rate, e.g. 0.1 = 10%: ";
	std::cin >> r;
	std::cout << "Option type, 1 = European call, 2 = Eurpoean put, 3 = Arithmetic Asian call, 4 = Arithmetic Asian put: ";
	std::cin >> type;

	OptionData myOption(k, T, r, sig, d, type);

	std::cout << "\nSimulation parameters\n";
	std::cout << "Start of underlying value: ";
	std::cin >> S_start;
	std::cout << "End of underlying value: ";
	std::cin >> S_end;
	std::cout << "grid size: ";
	std::cin >> grid_size;
	std::cout << "Number of simulations: ";
	std::cin >> NSIM;
	std::cout << "Number of time steps: ";
	std::cin >> NT;
	std::cout << "Type of simulation, Euler = 1 Exact1 = 2 Exact2 = 3: ";
	std::cin >> s;
	std::cout << std::endl;

	double ds = (S_end - S_start) / static_cast<double>(grid_size); // incremental increase in grid size
	std::vector<double> underlying(grid_size+1); // vector of underlying value

	// computing the underlying values
	underlying[0] = S_start;
	for (int i = 0; i < grid_size; i++)
	{
		underlying[i + 1] = underlying[i] + ds;
	}

	MonteCarlo MC(myOption, NSIM, NT, S_start);

	// ensuring the ouput on the console is clear
	const char separator = ' ';
	const int Width = 20;

	// Conditional to ensure the correct simulation type occurs i.e. s==1 euler, s==2 Exact1, s==3 Exact2, anytihing else returns "Incorrect choice of simulation type"
	if (s == 1)
	{
		std::vector<std::vector<double>> V = MC.Euler(S_start, S_end, grid_size);
		// Loop to ensure console output includes BS price if applicable i.e. BS price will be displayed if V[0].size() == 5
		if (V[0].size() == 5)
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "BS price" << "\n\n";

			for (int i = 0; i < grid_size+1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][4] << "\n";
			}
		}
		else
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time" << "\n\n";

			for (int i = 0; i < grid_size+1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3] << "\n";
			}
		}	
	}
	else if (s == 2)
	{
		std::vector<std::vector<double>> V = MC.Exact1(S_start, S_end, grid_size);
		if (V[0].size() == 5)
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "BS price" << "\n\n";

			for (int i = 0; i < grid_size + 1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][4] << "\n";
			}
		}
		else
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time" << "\n\n";

			for (int i = 0; i < grid_size + 1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3] << "\n";
			}
		}
	}
	else if (s == 3)
	{
		std::vector<std::vector<double>> V = MC.Exact2(S_start, S_end, grid_size);
		if (V[0].size() == 5)
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "BS price" << "\n\n";

			for (int i = 0; i < grid_size + 1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][4] << "\n";
			}
		}
		else
		{
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Underlying";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Simulated price";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Variance";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "Standard error";
			std::cout << std::left << std::setw(Width) << std::setfill(separator) << "time" << "\n\n";

			for (int i = 0; i < grid_size + 1; i++)
			{
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << underlying[i];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][0];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][1];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][2];
				std::cout << std::left << std::setw(Width) << std::setfill(separator) << V[i][3] << "\n";
			}
		}
	}
	else
	{
	std::cout << "Incorrect choice of simulation type" << std::endl;
	}

	return 0;
}