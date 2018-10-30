//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : CCDE.h
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Pawel Topa
//               Jaroslaw Was
// Version     : v1.0
// Created on  : Mar 20, 2016
//
// More details in the following paper:
//
// Trunfio, G.A., Topa, P., Was, J. 'A New Algorithm for Adapting the Configuration
// of Subcomponents in Large-Scale Optimization with Cooperative Coevolution, submitted'
//
//=======================================================================================

#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include "Benchmarks2010.h"
#include "JADE.h"
#include "Decomposer.h"
#include "numeric"

using namespace std;


class ConvPlotPoint
{
public:
    unsigned nfe;
    double f;
    unsigned subcomponentSize;
    unsigned individuals;
    ConvPlotPoint(unsigned  _nfe, double _f, unsigned _subcomponentSize, unsigned _individuals) :
        nfe(_nfe), f(_f), subcomponentSize(_subcomponentSize), individuals(_individuals)
    {};
};




/**
	@brief Main CCDE class.
	This class represents multiple swarms of particles which operate on the grouped directions of the main search space.
*/
class CCDE
{
    Fitness *fitness;

    ///Pseudorandom generator
    RandomEngine eng;

    uniform_real<double> unifRandom;

public:
    ///Create the CCDE object with the specified optimization parameters
    ///@param pNum number of particles
    ///@param D problem dimension
    ///@param ite number of generations
    ///@param
    CCDE();

    ///Destroy the CCDE object
    ~CCDE();

    ///Perform the optimization
    void optimize(Fitness* _function, unsigned int maxNumberOfEvaluations, unsigned sizeOfSubcomponents,
                  unsigned individualsPerSubcomponent, unsigned seed, unsigned nGenPerCycle);

    void optimize_MLSOFT(Fitness* _function, unsigned int maxNumberOfEvaluations, vector<unsigned> &_subcomponentSizes, unsigned numIndividualsPerSubcomponents,
                         unsigned seed, double tau, unsigned nGenPerCycle);


    void optimize_CCAS(Fitness* _function, unsigned int maxNumberOfEvaluations, vector<unsigned> &subcomponentSizes,
                       vector<unsigned> &numIndividualsPerSubcomponents, unsigned seed, unsigned nGenPerCycle, unsigned LCP,
                       double rho, double lambda, double delta, unsigned sigma);


    ///Returns the final fitness value
    double getFinalFitnessValue();

    double slope(const std::vector<double>& x, const std::vector<double>& y);

    ///Returns a pointer to the final global best position
    double* getFinalGlobalBestPosition();

    ///Prints results
    void printResults();

    void initPopulation(unsigned numOfIndividuals);

    void initContextVector();

    double computeFitness(vector<double> &x);

    Decomposer *createDecomposer(unsigned sizeOfSubcomponents, unsigned individualsPerSubcomponent, bool random = false);
    Decomposer *createDecomposer(unsigned indexOfDecomposerCharacteristics, bool random = false);

    void createSetOfDecomposers(vector<unsigned> &sizesOfSubcomponents, vector<unsigned> &individualsPerSubcomponent, bool random = false);

    unsigned findDecomposerWithMaxValueFunction();

    unsigned selectDecomposerSoftMax(double tau, vector<Decomposer*> &dec);

    unsigned selectDecomposerSoftMax(double tau, vector<double> &valueFunction);

    void copySearchState(unsigned bestFitnessDecomposer, unsigned selectedDecomposer);

    void selectNextPopulation(unsigned selectedDecomposer, unsigned bestFitnessDecomposer);

    void selectNextPopulation(unsigned selectedDecomposer);

    void broadcastSearchState(unsigned sourceDecomposer, double randomRatio);

    void optimizeSubcomponents(Decomposer *dec, unsigned nGenPerIteration);

    typedef enum { avgFitnessGain, stdDev, valueFunction } typeOfDecSort;
    struct doCompareDecomposers
    {
        doCompareDecomposers() { };

        bool operator()(const Decomposer* d1, const Decomposer * d2)
        {
            return d1->valueFunction > d2->valueFunction;
        }
    };

    ///Dimensionality of the search space
    unsigned problemDimension;

    //Set of decomposers
    vector< Decomposer* > decomposers;

    ///Number of fitness evaluations
    unsigned numberOfEvaluations;

    //JADE parameters
    double JADE_c;
    double JADE_p;
    int JADE_mutationStrategy;

    ///Final global best fitness value
    double globalBestFitness;

    ///Final global best position and context vector
    vector<double> contextVector;

    ///Current population
    vector< vector<double> > population;

    ///Fitnesses of population (only for initialization after sub-groups change)
    vector< double > fitnessValues;

    ///Lower limit for each dimension of the search space
    double lowerLimit;

    ///Upper limit for each dimension of the search space
    double upperLimit;

    ///index of the best individual
    //unsigned int globalBestIndex;

    ///Last elapsed time
    double elapsedTime;
    unsigned functionIndex;
    unsigned maxPopSize;
    vector<unsigned> subcomponentSizes;
    vector<unsigned> numIndividualsPerSubcomponents;
    unsigned maxNumberOfEvaluations;
};

