#include "MonteCarlo.hpp"
#include "OptionCommand.hpp"

int main()
{
	double S_0 = 25; // Spot price
	long NSIM = 1000; // Number of simulations
	long NT = 1000; // Number of time steps

	// imputting parameters of the option using boost library
	OptionData myOption((OptionParams::strike = 25, OptionParams::expiration = 0.25,
		OptionParams::volatility = 0.3, OptionParams::dividend = 0.0,
		OptionParams::optionType = 1, OptionParams::interestRate = 0.08));

	int M = 1000; // number of times running the member functions
	MonteCarlo Mc(myOption, NSIM, NT, S_0);
	std::vector<double> avgvalue_euler(4); //vector of the average value of the member function "Euler" output
	std::vector<double> avgvalue_exact1(4); //vector of the average value of the member function "Exact1" output
	std::vector<double> avgvalue_exact2(4); //vector of the average value of the member function "Exact2" output

	// caling the "euler" member function M times
	for (int i = 1; i <= M; i++)
	{
		std::vector<double> Value = Mc.Euler();
		avgvalue_euler[0] += Value[0];
		avgvalue_euler[1] += Value[1];
		avgvalue_euler[2] += Value[2];
		avgvalue_euler[3] += Value[3];
	}
	// getting the average of the "euler" member function ouput
	avgvalue_euler[0] = avgvalue_euler[0] / M;
	avgvalue_euler[1] = avgvalue_euler[1] / M;
	avgvalue_euler[2] = avgvalue_euler[2] / M;
	avgvalue_euler[3] = avgvalue_euler[3] / M;

	std::cout << "\nEuler:\nPrice: " << avgvalue_euler[0] << "\n" << "Standard Deviation: " << avgvalue_euler[1] << std::endl;
	std::cout << "Standard Error: " << avgvalue_euler[2] << "\n" << "Time: " << avgvalue_euler[3] << std::endl;

	// caling the "exact1" member function M times
	for (int i = 1; i <= M; i++)
	{
		std::vector<double> Value = Mc.Exact1();
		avgvalue_exact1[0] += Value[0];
		avgvalue_exact1[1] += Value[1];
		avgvalue_exact1[2] += Value[2];
		avgvalue_exact1[3] += Value[3];
	}
	// getting the average of the "exact1" member function ouput
	avgvalue_exact1[0] = avgvalue_exact1[0] / M;
	avgvalue_exact1[1] = avgvalue_exact1[1] / M;
	avgvalue_exact1[2] = avgvalue_exact1[2] / M;
	avgvalue_exact1[3] = avgvalue_exact1[3] / M;

	std::cout << "\nExact1:\nPrice: " << avgvalue_exact1[0] << "\n" << "Standard Deviation: " << avgvalue_exact1[1] << std::endl;
	std::cout << "Standard Error: " << avgvalue_exact1[2] << "\n" << "Time: " << avgvalue_exact1[3] << std::endl;

	// caling the "exact2" member function M times
	for (int i = 1; i <= M; i++)
	{
		std::vector<double> Value = Mc.Exact1();
		avgvalue_exact2[0] += Value[0];
		avgvalue_exact2[1] += Value[1];
		avgvalue_exact2[2] += Value[2];
		avgvalue_exact2[3] += Value[3];
	}
	// getting the average of the "exact1" member function ouput
	avgvalue_exact2[0] = avgvalue_exact2[0] / M;
	avgvalue_exact2[1] = avgvalue_exact2[1] / M;
	avgvalue_exact2[2] = avgvalue_exact2[2] / M;
	avgvalue_exact2[3] = avgvalue_exact2[3] / M;

	std::cout << "\nExact2:\nPrice: " << avgvalue_exact2[0] << "\n" << "Standard Deviation: " << avgvalue_exact2[1] << std::endl;
	std::cout << "Standard Error: " << avgvalue_exact2[2] << "\n" << "Time: " << avgvalue_exact2[3] << "\n" <<std::endl;

	// Displaying the Black Scholes call price
	CallPrice balckcall = CallPrice(myOption.K, myOption.T, myOption.r, myOption.D, myOption.sig);
	std::cout << "Black Scholes call Price: " << balckcall.execute(S_0) << std::endl;

	return 0;
}