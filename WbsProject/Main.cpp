#include "MonteCarlo.hpp"
#include "OptionCommand.hpp"

int main()
{
	double k;
	double T;
	double sig;
	double d;
	int type;
	double r;
	long NSIM;
	long NT;
	double S_0;
	bool a;

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
	std::cout << "Inital underlying value: ";
	std::cin >> S_0;
	std::cout << "Number of simulations: ";
	std::cin >> NSIM;
	std::cout << "Number of time steps: ";
	std::cin >> NT;
	std::cout << std::endl;
	
	// loop to produce Black Scholes price if applicable
	if (type == 2)
	{
		PutPrice balckput = PutPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
		std::cout << "Black Scholes put price: " << balckput.execute(S_0) << std::endl;
	}
	else if (type == 3)
	{
		std::cout << "No closed form solution\n";
	}
	else if (type == 4)
	{
		std::cout << "No closed form solution\n";
	}
	else 
	{
		CallPrice balckcall = CallPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
		std::cout << "Black Scholes call price: " << balckcall.execute(S_0) << std::endl;
	}

	std::cout << std::endl;

	// Outputting simulated values

	MonteCarlo Mc(myOption, NSIM, NT, S_0);
	std::vector<double> Value = Mc.Euler();
	std::cout << "Euler:\nPrice: " << Value[0] << "\n" << "Standard Deviation: " << Value[1] << std::endl;
	std::cout << "Standard Error: " << Value[2] << "\n" << "Time: " << Value[3] << std::endl;

	std::vector<double> Value1 = Mc.Exact1();
	std::cout << "\nExact1:\nPrice: " << Value1[0] << "\n" << "Standard Deviation: " << Value1[1] << std::endl;
	std::cout << "Standard Error: " << Value1[2] << "\n" << "Time: " << Value1[3] << std::endl;

	std::vector<double> Value2 = Mc.Exact2();
	std::cout << "\nExact2:\nPrice: " << Value2[0] << "\n" << "Standard Deviation: " << Value2[1] << std::endl;
	std::cout << "Standard Error: " << Value2[2] << "\n" << "Time: " << Value2[3] << std::endl;

	// generates new paths and saves these as an excel file if necissary
	std::cout << "\nSave an excel file gnerated by these finite difference methods, 1 = yes 0 = no: ";
	std::cin >> a;

	if (a == 1) 
	{
		Mc.Euler_excel();
		Mc.Exact1_excel();
		Mc.Exact2_excel();
		
		std::cout << "\nEnd of simulations\n";
	}
	else
	{
		std::cout << "\nEnd of simulations\n";
	}


	return 0;

}

