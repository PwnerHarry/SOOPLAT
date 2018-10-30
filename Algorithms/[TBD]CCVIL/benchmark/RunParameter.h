#ifndef _RUNPARAMETER_H
#define _RUNPARAMETER_H

#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include	<cmath>
#include	<ctime>

using namespace std;

class RunParameter{
	public:
		// Dimension of problems
		unsigned dimension;

		// the amount of independent run
		unsigned numOfRun;

		// initial population size
		unsigned NP;

		// initial Group Size
		unsigned initialGroupSize;

		// The number of Sampling points for plotting the convergence curve
		unsigned samplingPoint;

		// Sampling interval for plotting the convergence curve
		unsigned samplingInterval;

		// initialized random seed
		int initRandomSeed;

		// group size for non-separable part of function
		unsigned nonSeparableGroupSize;

		// strategy for learning and grouping 
		unsigned learnStrategy;

		// Fitness check point
		vector<unsigned> fitnessCheckPoint;

		// the IDes of benchmark functions to be tested in the experiment
		vector<unsigned> functionToRun;

		// CCVIL's lower threshold
		unsigned lowerThreshold;

		// for running the partial information experiment
		// the percentage of prior knonw grouping information
		vector<double> knownGroupPercent;

		double learnPortion;

		// Runtime Parameter for rJADE
		double c;	

		double p;

		unsigned failThreshold;

		unsigned Afactor;

		unsigned performOpt; 

		// default constructor
		RunParameter();

		// constructor with configuration file name specified
		RunParameter(char* filePath); 

		// default destructor
		~RunParameter();
};
#endif
