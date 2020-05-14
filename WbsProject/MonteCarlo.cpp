// cpp file for class of montecarlo

#include "MonteCarlo.hpp"

// Constructors & Destructors 
MonteCarlo::MonteCarlo() : myOption(new OptionData(0,0,0,0,0,1)), NSIM(0), NT(0), S_0(0) {} // deafult constructor
MonteCarlo::MonteCarlo(const OptionData& optionData, long Sim, long Time, double Spot) :
    myOption(new OptionData(optionData)), NSIM(Sim), NT(Time), S_0(Spot) {};
MonteCarlo::MonteCarlo(const MonteCarlo& Mc) : myOption(Mc.myOption), NSIM(Mc.NSIM), NT(Mc.NT), S_0(Mc.S_0) {}; // copy constructor
MonteCarlo::~MonteCarlo() {}; //destructor

// Member functions 

// Euler simulations
std::vector<double> MonteCarlo::Euler()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);
	
	double price = 0.0;	// Option price
	double payoffT; // payoff of option at T
	double squaredPayoff = 0.0; // squared payoff of option at T
	double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	//starting the timer
	StopWatch<> sw;
	sw.Start();
	// preconditioning all the random variates for each path and each simulation
	for (int i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for
	for (long i = 1; i <= NSIM; i++)
	{ 
		std::vector<double> path(NT + 1); // Creating vector of timesteps

		path[0] = S_0;  // Setting the first element of the path equal to S_0
		x = 0.0; // Setting the time elapsed to 0
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case explicit Euler)
			path[index + 1] = path[index] + (dt * sde.drift_euler(x, path[index])) + (sqrdt * sde.diffusion_euler(x, path[index]) * dW[i-1][index]);

			// increasing time elapsed
			x += dt;
		}

		paths[i - 1] = path; // updating the vector of simulations

		// Assemble quantities 
		payoffT = myOption->myPayOffFunction(paths[i - 1]);
		sumPriceT += payoffT;
		squaredPayoff += (payoffT * payoffT);
	}
	// Finally, discounting the average payoff
	price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

	// Computing the standard deviation and error
	double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
	double SE = SD / std::sqrt(M);

	//stopping the stopwatch and getting the time elapsed
	sw.Stop();
	double time = sw.GetTime();

	// creating and populating the outputted vector of price, standard deviation, standard error and time elapsed
	std::vector<double> V;
	V.push_back(price);
	V.push_back(SD);
	V.push_back(SE);
	V.push_back(time);

	// loop to produce Black Scholes price if applicable
	if (myOption->type == 2)
	{
		V.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}
	else if (myOption->type == 3)
	{

	}
	else if (myOption->type == 4)
	{

	}
	else
	{
		V.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}

	return V;

};

// Exact Simulation version 1
std::vector<double> MonteCarlo::Exact1()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	double price = 0.0;	// Option price
	double payoffT; // payoff of option at T
	double squaredPayoff = 0.0; // squared payoff of option at T
	double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	// starting the timer
	StopWatch<> sw;
	sw.Start();
	// preconditioning all the random variates for each path and each simulation
	for (long i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for
	for (long i = 1; i <= NSIM; i++)
	{

		std::vector<double> path(NT + 1);  // Creating vector of timesteps
		x = 0.0; // Setting the time elapsed to 0
		path[0] = S_0; // Setting the first element of the path equal to S_0
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case exact simulation verion 1)
			path[index + 1] = path[index] * std::exp(dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i-1][index]);

			// increasing time elapsed
			x += dt;

		}

		paths[i - 1] = path; // updating the vector of simulations

		// Assemble quantities
		payoffT = myOption->myPayOffFunction(paths[i - 1]);
		sumPriceT += payoffT;
		squaredPayoff += (payoffT * payoffT);

	}
	// Finally, discounting the average price
	price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

	// Computing the standard deviation and error
	double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
	double SE = SD / std::sqrt(M);

	//stopping the stopwatch and getting the time elapsed
	sw.Stop();
	double time = sw.GetTime();

	// creating and populating the outputted vector of price, standard deviation, standard error and time elapsed
	std::vector<double> V;
	V.push_back(price);
	V.push_back(SD);
	V.push_back(SE);
	V.push_back(time);
	
	if (myOption->type == 2)
	{
		V.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}
	else if (myOption->type == 3)
	{

	}
	else if (myOption->type == 4)
	{

	}
	else
	{
		V.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}

	return V;


}

// Exact Simulation version 2
std::vector<double> MonteCarlo::Exact2()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	double price = 0.0;	// Option price
	double payoffT; // payoff of option at T
	double squaredPayoff = 0.0; // squared payoff of option at T
	double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	// starting the timer
	StopWatch<> sw;
	sw.Start();
	// preconditioning all the random variates for each path and each simulation
	for (long i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for
	for (long i = 1; i <= NSIM; i++)
	{

		std::vector<double> path(NT + 1); // Creating vector of timesteps
		std::vector<double> coefficent(NT + 1); // Creating vector of "coeffiecents" given as x in the project
		x = 0.0; // Setting the time elapsed to 0
		path[0] = S_0; // Setting the first element of the path equal to S_0
		coefficent[0] = std::log(S_0); // Setting the first value of the coefficent vector
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case explicit Exact simulation version 2)
			coefficent[index + 1] = coefficent[index] + dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i-1][index];
			path[index + 1] = std::exp(coefficent[index]);

			// increasing time elapsed
			x += dt;
		}

		paths[i - 1] = path; // updating the vector of simulations

		// Assemble quantities
		payoffT = myOption->myPayOffFunction(paths[i - 1]);
		sumPriceT += payoffT;
		squaredPayoff += (payoffT * payoffT);
	}
	// Finally, discounting the average price
	price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

	// Computing the standard deviation and error
	double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
	double SE = SD / std::sqrt(M);

	//stopping the stopwatch and getting the time elapsed
	sw.Stop();
	double time = sw.GetTime();

	// creating and populating the outputted vector of price, standard deviation, standard error and time elapsed
	std::vector<double> V;
	V.push_back(price);
	V.push_back(SD);
	V.push_back(SE);
	V.push_back(time);

	if (myOption->type == 2)
	{
		V.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}
	else if (myOption->type == 3)
	{

	}
	else if (myOption->type == 4)
	{

	}
	else
	{
		V.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S_0));
	}

	return V;

}

void MonteCarlo::Euler_excel()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	// preconditioning all the random variates for each path and each simulation
	for (long i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for 
	for (long i = 1; i <= NSIM; i++)
	{

		std::vector<double> path(NT + 1); // Creating vector of timesteps

		path[0] = S_0;  // Setting the first element of the path equal to S_0
		x = 0.0; // Setting the time elapsed to 0
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case explicit Exact simulation version 2)
			path[index + 1] = path[index] + (dt * sde.drift_euler(x, path[index])) + (sqrdt * sde.diffusion_euler(x, path[index]) * dW[i-1][index]);

			// increasing time elapsed
			x += dt;
		}

		paths[i - 1] = path; // updating the vector of simulations

	}
	//produce csv file
	data_output_euler(paths, NSIM, NT);

};

void MonteCarlo::Exact1_excel()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	// preconditioning all the random variates for each path and each simulation
	for (long i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for 
	for (long i = 1; i <= NSIM; i++)
	{ // Calculate a path at each iteration
		std::vector<double> path(NT + 1); // Creating vector of timesteps

		path[0] = S_0;  // Setting the first element of the path equal to S_0
		x = 0.0; // Setting the time elapsed to 0
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case explicit Exact simulation version 1)
			path[index + 1] = path[index] * std::exp(dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i-1][index]);

			// increasing time elapsed
			x += dt;

		}

		paths[i - 1] = path; // updating the vector of simulations

	}

	//produce csv file
	data_output_exact1(paths, NSIM, NT);

}

void MonteCarlo::Exact2_excel()
{
	//generating SDE class from Option data 
	SDE sde(*myOption.get());

	// x created to keep track of the time elapsed
	double x;

	// computing size of indiviudal time step and sqrt of time step
	double M = static_cast<double>(NSIM);
	double dt = myOption->T / double(NT);
	double sqrdt = std::sqrt(dt);

	// Vectors used to pre condition random variates
	std::vector<std::vector<double>> dW(NSIM);
	std::vector<double> m(NT + 1);

	std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

	// preconditioning all the random variates for each path and each simulation
	for (long i = 1; i <= M; i++)
	{
		generateRN(m);
		dW[i - 1] = m;
	}

	/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

#pragma omp parallel for 
	for (long i = 1; i <= NSIM; i++)
	{ // Calculate a path at each iteration
		std::vector<double> path(NT + 1); // Creating vector of timesteps
		std::vector<double> coefficent(NT + 1); // Creating vector of "coeffiecents" given as x in the project
		path[0] = S_0;  // Setting the first element of the path equal to S_0
		x = 0.0; // Setting the time elapsed to 0
		coefficent[0] = std::log(S_0);// Setting the first value of the coefficent vector
		for (long index = 0; index < NT; index++)
		{

			// The FDM (in this case exact simulation version 2)
			coefficent[index + 1] = coefficent[index] + dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i-1][index];
			path[index + 1] = std::exp(coefficent[index]);

			// increasing time elapsed
			x += dt;

		}

		paths[i - 1] = path; // updating the vector of simulations

	}

	//produce csv file
	data_output_exact2(paths, NSIM, NT);

}

//////// Member functions to produce prices for range of spots ////////

std::vector<std::vector<double>> MonteCarlo::Euler(double S_start, double S_end, int grid_size)
{
	// precondition checks for grid size
	if (grid_size < 1)
	{
		std::cout << "\nGrid size cannot be smaller than 1" << std::endl;
	}
	else
	{
		// Creating vector of spot price values
		std::vector<double> S(grid_size + 1);
		S[0] = S_start;
		double g_size = static_cast<double>(grid_size);

		for (int i = 0; i < grid_size; i++)
		{
			S[i + 1] = S[i] + (S_end - S_start) / g_size;
		}

		// Output vector which will have Price, Variance, Standard Error, and time for all spot prices
		std::vector<std::vector<double>> V;

		//generating SDE class from Option data 
		SDE sde(*myOption.get());

		// computing size of indiviudal time step and sqrt of time step
		double M = static_cast<double>(NSIM);
		double dt = myOption->T / double(NT);
		double sqrdt = std::sqrt(dt);

		// Vectors used to pre condition random variates
		std::vector<std::vector<double>> dW(NSIM);
		std::vector<double> m(NT + 1);

		std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

		// preconditioning all the random variates for each path and each simulation
		for (long i = 1; i <= M; i++)
		{
			generateRN(m);
			dW[i - 1] = m;
		}

		/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

		for (int j = 0; j <= grid_size; j++)
		{
			// x created to keep track of the time elapsed
			double x;

			double price = 0.0;	// Option price
			double payoffT; // payoff of option at T
			double squaredPayoff = 0.0; // squared payoff of option at T
			double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

			//starting the timer
			StopWatch<> sw;
			sw.Start();
#pragma omp parallel for
			for (long i = 1; i <= NSIM; i++)
			{
				std::vector<double> path(NT + 1); // Creating vector of timesteps

				path[0] = S[j];  // Setting the first element of the path equal to S[j]
				x = 0.0; // Setting the time elapsed to 0
				for (long index = 0; index < NT; index++)
				{

					// The FDM (in this case explicit Euler)
					path[index + 1] = path[index] + (dt * sde.drift_euler(x, path[index])) + (sqrdt * sde.diffusion_euler(x, path[index]) * dW[i - 1][index]);

					// increasing time elapsed
					x += dt;
				}

				paths[i - 1] = path; // updating the vector of simulations

				// Assemble quantities 
				payoffT = myOption->myPayOffFunction(paths[i - 1]);
				sumPriceT += payoffT;
				squaredPayoff += (payoffT * payoffT);
			}
			// Finally, discounting the average payoff
			price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

			// Computing the standard deviation and error
			double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
			double SE = SD / std::sqrt(M);

			//stopping the stopwatch and getting the time elapsed
			sw.Stop();
			double time = sw.GetTime();

			std::vector<double> v_element;
			v_element.push_back(price);
			v_element.push_back(SD);
			v_element.push_back(SE);
			v_element.push_back(time);
			
			// loop to produce Black Scholes price if applicable
			if (myOption->type == 2)
			{
				v_element.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}
			else if (myOption->type == 3)
			{
			
			}
			else if (myOption->type == 4)
			{
			
			}
			else
			{
				v_element.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}

			V.push_back(v_element);
		}

		// Loop to ensure the excel file contains a closed form solution to the option price if applicable
		if (myOption->type == 2)
		{
			data_output_euler_range(V, S, 5);
		}
		else if (myOption->type == 3)
		{
			data_output_euler_range(V, S, 4);
		}
		else if (myOption->type == 4)
		{
			data_output_euler_range(V, S, 4);
		}
		else
		{
			data_output_euler_range(V, S, 5);;
		}

		return V;
	}

}

std::vector<std::vector<double>> MonteCarlo::Exact1(double S_start, double S_end, int grid_size)
{
	// precondition checks for grid size
	if (grid_size < 1)
	{
		std::cout << "\nGrid size cannot be smaller than 1" << std::endl;
	}
	else
	{
		// Creating vector of spot price values
		std::vector<double> S(grid_size + 1);
		S[0] = S_start;
		double g_size = static_cast<double>(grid_size);

		for (int i = 0; i < grid_size; i++)
		{
			S[i + 1] = S[i] + (S_end - S_start) / g_size;
		}

		// Output vector which will have Price, Variance, Standard Error, and time for all spot prices
		std::vector<std::vector<double>> V;

		//generating SDE class from Option data 
		SDE sde(*myOption.get());

		// computing size of indiviudal time step and sqrt of time step
		double M = static_cast<double>(NSIM);
		double dt = myOption->T / double(NT);
		double sqrdt = std::sqrt(dt);

		// Vectors used to pre condition random variates
		std::vector<std::vector<double>> dW(NSIM);
		std::vector<double> m(NT + 1);

		std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

		// preconditioning all the random variates for each path and each simulation
		for (long i = 1; i <= M; i++)
		{
			generateRN(m);
			dW[i - 1] = m;
		}

		/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

		for (int j = 0; j <= grid_size; j++)
		{
			// x created to keep track of the time elapsed
			double x;

			double price = 0.0;	// Option price
			double payoffT; // payoff of option at T
			double squaredPayoff = 0.0; // squared payoff of option at T
			double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

			//starting the timer
			StopWatch<> sw;
			sw.Start();
#pragma omp parallel for
			for (long i = 1; i <= NSIM; i++)
			{
				std::vector<double> path(NT + 1); // Creating vector of timesteps

				path[0] = S[j];  // Setting the first element of the path equal to S[j]
				x = 0.0; // Setting the time elapsed to 0
				for (long index = 0; index < NT; index++)
				{

					// The FDM (in this case exact simulation verion 1)
					path[index + 1] = path[index] * std::exp(dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i - 1][index]);

					// increasing time elapsed
					x += dt;
				}

				paths[i - 1] = path; // updating the vector of simulations

				// Assemble quantities 
				payoffT = myOption->myPayOffFunction(paths[i - 1]);
				sumPriceT += payoffT;
				squaredPayoff += (payoffT * payoffT);
			}
			// Finally, discounting the average payoff
			price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

			// Computing the standard deviation and error
			double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
			double SE = SD / std::sqrt(M);

			//stopping the stopwatch and getting the time elapsed
			sw.Stop();
			double time = sw.GetTime();

			std::vector<double> v_element;
			v_element.push_back(price);
			v_element.push_back(SD);
			v_element.push_back(SE);
			v_element.push_back(time);

			// loop to produce Black Scholes price if applicable
			if (myOption->type == 2)
			{
				v_element.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}
			else if (myOption->type == 3)
			{

			}
			else if (myOption->type == 4)
			{

			}
			else
			{
				v_element.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}

			V.push_back(v_element);
		}

		// Loop to ensure the excel file contains a closed form solution to the option price if applicable
		if (myOption->type == 2)
		{
			data_output_exact1_range(V, S, 5);
		}
		else if (myOption->type == 3)
		{
			data_output_exact1_range(V, S, 4);
		}
		else if (myOption->type == 4)
		{
			data_output_exact1_range(V, S, 4);
		}
		else
		{
			data_output_exact1_range(V, S, 5);;
		}

		return V;
	}
}


std::vector<std::vector<double>> MonteCarlo::Exact2(double S_start, double S_end, int grid_size)
{
	// precondition checks for grid size
	if (grid_size < 1)
	{
		std::cout << "\nGrid size cannot be smaller than 1" << std::endl;
	}
	else
	{
		// Creating vector of spot price values
		std::vector<double> S(grid_size + 1);
		S[0] = S_start;
		double g_size = static_cast<double>(grid_size);

		for (int i = 0; i < grid_size; i++)
		{
			S[i + 1] = S[i] + (S_end - S_start) / g_size;
		}

		// Output vector which will have Price, Variance, Standard Error, and time for all spot prices
		std::vector<std::vector<double>> V;

		//generating SDE class from Option data 
		SDE sde(*myOption.get());

		// computing size of indiviudal time step and sqrt of time step
		double M = static_cast<double>(NSIM);
		double dt = myOption->T / double(NT);
		double sqrdt = std::sqrt(dt);

		// Vectors used to pre condition random variates
		std::vector<std::vector<double>> dW(NSIM);
		std::vector<double> m(NT + 1);

		std::vector<std::vector<double>> paths(NSIM); // Creating 2D vector of simulations

		// preconditioning all the random variates for each path and each simulation
		for (long i = 1; i <= M; i++)
		{
			generateRN(m);
			dW[i - 1] = m;
		}

		/* // Set the number of threads to use (optional)
		omp_set_num_threads(4);	// Number of threads = 4 */

		for (int j = 0; j <= grid_size; j++)
		{
			// x created to keep track of the time elapsed
			double x;

			double price = 0.0;	// Option price
			double payoffT; // payoff of option at T
			double squaredPayoff = 0.0; // squared payoff of option at T
			double sumPriceT = 0.0; // sum of payoffs at time T for all simulations

			//starting the timer
			StopWatch<> sw;
			sw.Start();
#pragma omp parallel for
			for (long i = 1; i <= NSIM; i++)
			{
				std::vector<double> path(NT + 1); // Creating vector of timesteps
				std::vector<double> coefficent(NT + 1); // Creating vector of "coefficents" given as x in the project

				path[0] = S[j];  // Setting the first element of the path equal to S[j]
				x = 0.0; // Setting the time elapsed to 0
				coefficent[0] = std::log(S[j]); // Setting the first value of the coefficent vector
				for (long index = 0; index < NT; index++)
				{

					// The FDM (in this case explicit Exact simulation version 2)
					coefficent[index + 1] = coefficent[index] + dt * sde.deterministic_exact() + myOption->sig * sqrdt * dW[i - 1][index];
					path[index + 1] = std::exp(coefficent[index]);

					// increasing time elapsed
					x += dt;
				}

				paths[i - 1] = path; // updating the vector of simulations

				// Assemble quantities 
				payoffT = myOption->myPayOffFunction(paths[i - 1]);
				sumPriceT += payoffT;
				squaredPayoff += (payoffT * payoffT);
			}
			// Finally, discounting the average payoff
			price = std::exp(-myOption->r * myOption->T) * sumPriceT / M;

			// Computing the standard deviation and error
			double SD = std::sqrt((squaredPayoff / M) - sumPriceT * sumPriceT / (M * M));
			double SE = SD / std::sqrt(M);

			//stopping the stopwatch and getting the time elapsed
			sw.Stop();
			double time = sw.GetTime();

			std::vector<double> v_element;
			v_element.push_back(price);
			v_element.push_back(SD);
			v_element.push_back(SE);
			v_element.push_back(time);

			// loop to produce Black Scholes price if applicable
			if (myOption->type == 2)
			{
				v_element.push_back(PutPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}
			else if (myOption->type == 3)
			{

			}
			else if (myOption->type == 4)
			{

			}
			else
			{
				v_element.push_back(CallPrice(myOption->K, myOption->T, myOption->r, myOption->D, myOption->sig).execute(S[j]));
			}

			V.push_back(v_element);
		}

		// Loop to ensure the excel file contains a closed form solution to the option price if applicable
		if (myOption->type == 2)
		{
			data_output_exact2_range(V, S, 5);
		}
		else if (myOption->type == 3)
		{
			data_output_exact2_range(V, S, 4);
		}
		else if (myOption->type == 4)
		{
			data_output_exact2_range(V, S, 4);
		}
		else
		{
			data_output_exact2_range(V, S, 5);;
		}

		return V;
	}
}
