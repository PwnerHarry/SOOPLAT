//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : Decomposer.h
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
#include <vector>
#include <map>
#include <random>
#include <functional>
#include <numeric>
#include "JADE.h"

class CCDE;

using namespace std;

typedef mt19937 RandomEngine;

class Decomposer
{
public:
    typedef enum { fromBestToWorst, fromWorstToBest} typeOfSort;

    struct Reward
    {
        Reward(double _r, unsigned _iteration) : r(_r), iteration(_iteration) {};
        double r;
        unsigned iteration;
    };

    vector<unsigned> coordinates;
    map<unsigned, JADE*> coordinateToOptimizer;
    CCDE &CCOptimizer;
    vector< JADE* > optimizers;
    int sizeOfSubcomponents;
    int indexOfDecomposerCharacteristics;
    bool applyRandomGrouping;
    double bestAchievedFitness;
    double valueFunction;
    vector <Reward> rewards;
    RandomEngine eng;
    unsigned individualsPerSubcomponent;
    double prevBestFitness;
    unsigned functionEvaluations;
    unsigned learningCost;
    double prob;
    unsigned costOfCycle;
    unsigned learningCycles;

    //Current population
    vector< vector<double> > population;

    //Final global best position and context vector
    vector<double> contextVector;

    //Fitnesses of population
    vector< double > fitnessValues;

    vector< unsigned> popSort;

    struct doCompareDecomposers
    {
        vector<Decomposer*> &groups;
        bool sortByFitness;
        doCompareDecomposers(vector<Decomposer*> &g, bool _sortByFitness = false) :
            groups(g), sortByFitness(_sortByFitness) { };

        bool operator()(const unsigned g1i, const unsigned g2i)
        {
            Decomposer *g1 = groups[g1i];
            Decomposer *g2 = groups[g2i];
            if ( !sortByFitness )
            {
                //dal migliore al peggiore
                return g1->valueFunction > g2->valueFunction;
            }
            else
            {
                //dal peggiore al migliore
                return g1->bestAchievedFitness > g2->bestAchievedFitness;
            }
        }
    };

    struct doCompareIndividuals
    {
        vector<double> &fitness;
        typeOfSort type;
        doCompareIndividuals(vector<double> &_fitness, typeOfSort _type) : fitness(_fitness), type(_type) { };

        bool operator()(const unsigned i1, const unsigned i2)
        {
            if ( type==fromBestToWorst )
                return fitness[i1] < fitness[i2];
            else
                return fitness[i1] > fitness[i2];
        }
    };

    Decomposer(CCDE &_CCOptimizer, unsigned seed, vector<unsigned> &_coordinates,
               unsigned _sizeOfSubcomponents,
               unsigned _individualsPerSubcomponent,
               vector< vector<double> > &_population,
               vector<double>  &_contextVector,
               bool RG);
    ~Decomposer();
    void setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, vector<double> &fitnessValues);
    void setCoordinates(vector<unsigned> &_coordinates);
    void updateContextVector(JADE *optimizer);
    void buildContextVector();
    void randomGrouping();
    void sortPopulation(typeOfSort t);
    void setSeed(unsigned seed);
    void setOptimizersCoordinatesAndEvaluatePopulation();


};
