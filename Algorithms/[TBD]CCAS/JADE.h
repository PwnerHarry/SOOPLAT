//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : JADE.h
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
#include <map>
#include <cstddef>


class Decomposer;

using namespace std;

///To chose between maximisation and minimisation
//#define MAXIMIZE
#ifdef MAXIMIZE
#define BETTER_THAN >
#else
#define BETTER_THAN <
#endif

class JADE
{
    Decomposer &decomposer;

    struct doCompareIndividuals
    {
        doCompareIndividuals(const double *_f) : f(_f) { }
        const double *f;

        bool operator()(const int & i1, const int & i2)
        {
            return f[i1] < f[i2];
        }
    };

public:
    JADE(unsigned _dimension, unsigned _numberOfIndividuals, Decomposer &_group);
    void setCoordinates(vector<unsigned> &_coordinates);
    void setCoordinates(unsigned *coordinates, unsigned numOfCoordinates);

    void update();
    void sortPopulation(vector<double> &fitness, vector<int> &sortIndex);
    void evaluatePopulation(vector< vector<double> > &population, vector< double > &fitness);
    int evaluateParents();
    double calculateFitnessValue(vector<double> &p);
    int optimize(int iterations);
    void updateIndexOfBest();
    void loadIndividuals(vector< vector<double> > &population);
    void storeIndividuals(vector< vector<double> > &population);
    void setParentFitness(vector<double> &fitnessValues);

    vector<double> &getCollaborator();
    vector<unsigned> coordinates;
    map<unsigned, unsigned> globalCoordToLocalCoord;
    unsigned int dimension;

    ///array containing the positions of all individuals
    vector< vector<double> > parents;
    vector< vector<double> > offsprings;

    vector<int> sortIndex;
    vector<double> FF;
    vector<double> CR;
    vector<double> SSFF;  // successful F values
    vector<double> SSCR;  // successful CR values
    vector<int> binaryVector;

    double	JADE_mu_cr;
    double	JADE_mu_ff;

    ///array containing the current fitness of all particles
    vector< double > parentsFitness;
    vector< double > offspringsFitness;

    double bestFitness;

    ///array containing the index of the best position attained so far
    unsigned indexOfBest;

    unsigned numberOfIndividuals;

    ///buffer
    vector< double > xp;

    uniform_real<double> unifRandom;
};

