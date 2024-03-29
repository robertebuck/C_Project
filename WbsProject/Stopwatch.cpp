// Taken directly from StopWatch.cpp given in class used as a black box no editing made
// with the exception of including #ifndef, #define, #endif to prevent multiple inclusions 

#ifndef Stopwatch_CPP
#define Stopwatch_CPP

#include"Stopwatch.hpp"
#include<iostream>

template <typename TickType,typename UnitType>
	void StopWatch<TickType,UnitType>::Start()
{
	start = std::chrono::system_clock::now();
	isStart = true;
}

template <typename TickType,typename UnitType>
	void StopWatch<TickType,UnitType>::Stop() 
{
	end = std::chrono::system_clock::now();
	isEnd = true;
}

template <typename TickType,typename UnitType>
	void StopWatch<TickType,UnitType>::Reset() 
{
	start = std::chrono::system_clock::now();
	isStart = false;
	end = std::chrono::system_clock::now();
	isEnd = false;
}

template <typename TickType,typename UnitType>
	TickType StopWatch<TickType,UnitType>::GetTime()  
{
	TickType count;
	if (isStart && isEnd) 
	{
		// Seconds as unit type (default)
		std::chrono::duration<TickType, UnitType> D = end - start;
		count = D.count(); 
	}
	else 
	{
		std::cout << "Start time point or end time point not recorded\n";
		count = -1;
	}

	// Reset all counters
	Reset();

	return count;
}

#endif