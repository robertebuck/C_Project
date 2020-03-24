#include "MonteCarlo.hpp"
#include "OptionCommand.hpp"

int main()
{
	double S_start;
	double S_end;
	double increments;
	double k;
	double T;
	double sig;
	double d;
	int type;
	double r;
	long NSIM;
	long NT;
	int s;
	std::vector<double> prices; // vector of prices for each spot
	std::vector<double> underlying; // vector of underlying value
	std::vector<double> bsprice; // vector of Black Scholes price for each underlying value

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
	std::cout << "Incremental increase in underlying value: ";
	std::cin >> increments;
	std::cout << "Number of simulations: ";
	std::cin >> NSIM;
	std::cout << "Number of time steps: ";
	std::cin >> NT;
	std::cout << "Type of simulation, Euler = 1 Exact1 = 2 Exact2 = 3: ";
	std::cin >> s;
	std::cout << std::endl;

	// for European put option i.e type == 2
	if (type == 2)
	{ // for "Exact1" i.e. s == 2 simulation
		if (s == 2)
		{ // looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				PutPrice blackput = PutPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Exact1()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackput.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact1.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}
		// for "Exact2" i.e. s == 3 simulation
		else if (s == 3)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				PutPrice blackput = PutPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Exact2()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackput.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact2.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}
		// for "Euler" simulation i.e. s == 1 or just s not = 2 or 3
		else
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				PutPrice blackput = PutPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Euler()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackput.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_euler.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}
	}
	// for Antithetic Asian call option i.e. type == 3
	else if (type == 3)
	{// for "Exact1" i.e. s == 2 simulation
		if (s == 2)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Exact1()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact1.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
		// for "Exact2" i.e. s == 3 simulation
		else if (s == 3)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Exact2()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact2.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
		// for "Euler" simulation i.e. s == 1 or just s not = 2 or 3
		else
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Euler()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_euler.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
	}
	// For Antithetic Asian put options i.e. type == 4
	else if (type == 4)
	{// for "Exact1" i.e. s == 2 simulation
		if (s == 2)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Exact1()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact1.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
		// for "Exact2" i.e. s == 3 simulation
		else if (s == 3)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Exact2()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact2.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
		// for "Euler" simulation i.e. s == 1 or just s not = 2 or 3
		else
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				prices.push_back(Mc.Euler()[0]);
				underlying.push_back(i);
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_euler.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << '\n';
			}
		}
	}
	// For European call option i.e. type == 1 or just type not == 2,3 or 4
	else
	{// for "Exact1" i.e. s == 2 simulation
		if (s == 2)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				CallPrice blackcall = CallPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Exact1()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackcall.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact1.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}
		// for "Exact2" i.e. s == 3 simulation
		else if (s == 3)
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				CallPrice blackcall = CallPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Exact2()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackcall.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_exact2.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}
		// for "Euler" simulation i.e. s == 1 or just s not = 2 or 3
		else
		{// looping through all the underlying values
			for (double i = S_start; i <= S_end; i += increments)
			{
				MonteCarlo Mc(myOption, NSIM, NT, i);
				CallPrice blackcall = CallPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
				prices.push_back(Mc.Euler()[0]);
				underlying.push_back(i);
				bsprice.push_back(blackcall.execute(i));
			}
			// saving output of underlying, simulated prices and the Black-Scholes price in an excel file
			std::ofstream out("Price_range_euler.csv");
			for (std::size_t j = 0; j < prices.size(); j++)
			{
				out << underlying[j] << ',' << prices[j] << ',' << bsprice[j] << '\n';
			}
		}

	}
	return 0;
}