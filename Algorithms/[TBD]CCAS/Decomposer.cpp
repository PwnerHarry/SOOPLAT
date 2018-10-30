//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : Decomposer.cpp
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


#include "Decomposer.h"
#include "CCDE.h"


Decomposer::Decomposer(CCDE &_CCOptimizer,
                       unsigned seed,
                       vector<unsigned> &_coordinates,
                       unsigned _sizeOfSubcomponents,
                       unsigned _individualsPerSubcomponent,
                       vector< vector<double> > &_population,
                       vector<double>  &_contextVector,
                       bool RG) :  CCOptimizer(_CCOptimizer), sizeOfSubcomponents(_sizeOfSubcomponents),
    individualsPerSubcomponent(_individualsPerSubcomponent), applyRandomGrouping(RG)
{
    eng.seed(seed);

    bestAchievedFitness = 1.0E100;
    coordinates = _coordinates;
    indexOfDecomposerCharacteristics = -1;

    for (unsigned i = 0; i < individualsPerSubcomponent; ++i)
        population.push_back(_population[i]);

    contextVector.resize(CCOptimizer.problemDimension);
    for (unsigned i = 0; i < CCOptimizer.problemDimension; ++i)
        contextVector[i] = _contextVector[i];

    unsigned d = 0, size = sizeOfSubcomponents;
    while ( d<coordinates.size() )
    {
        if ( d + size > coordinates.size() )
            size = coordinates.size() - d;

        JADE *optimizer = new JADE(size, individualsPerSubcomponent, *this);

        optimizer->setCoordinates(&(coordinates[d]), size);

        for (int k = 0; k < size; ++k)
            coordinateToOptimizer[coordinates[d + k]] = optimizer;

        optimizer->loadIndividuals(population);
        optimizer->evaluateParents();

        optimizers.push_back(optimizer);

        d += size;
    }

    popSort.clear();
    for (unsigned i = 0; i < population.size(); ++i)
        popSort.push_back(i);

}



void Decomposer::setOptimizersCoordinatesAndEvaluatePopulation()
{
    unsigned d = 0, k = 0, size = sizeOfSubcomponents;
    while ( d<coordinates.size() )
    {
        if (d + size > coordinates.size())
            size = coordinates.size() - d;

        JADE *optimizer = optimizers[k];

        optimizer->setCoordinates(&(coordinates[d]), size);

        for (int q = 0; q < size; ++q)
            coordinateToOptimizer[coordinates[d + q]] = optimizer;

        optimizer->loadIndividuals(population);
        optimizer->evaluateParents();
        optimizer->updateIndexOfBest();

        d += size;
        k++;
    }

}

Decomposer::~Decomposer()
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
};



void Decomposer::setSeed(unsigned seed)
{
    eng.seed(seed);
}



void Decomposer::setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, vector<double> &fitnessValues)
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
    optimizers.clear();

    sizeOfSubcomponents = newSizeOfSubcomponents;
    unsigned numberOfSubcomponents = coordinates.size() / newSizeOfSubcomponents;

    individualsPerSubcomponent = fitnessValues.size();

    for (unsigned i = 0; i<numberOfSubcomponents; ++i)
    {
        optimizers.push_back(new JADE(sizeOfSubcomponents, individualsPerSubcomponent, *this));

        optimizers[i]->setCoordinates(&(coordinates[i*sizeOfSubcomponents]), sizeOfSubcomponents);

        for (int k = 0; k < sizeOfSubcomponents; ++k)
            coordinateToOptimizer[coordinates[i*sizeOfSubcomponents + k]] = optimizers[i];

        optimizers[i]->loadIndividuals(population);
        optimizers[i]->setParentFitness(fitnessValues);

        //optimizers[i]->evaluateParents();

        // optimizers[i]->updateIndexOfBest();
    }
}



void Decomposer::setCoordinates(vector<unsigned> &_coordinates)
{
    coordinates = coordinates;
}



void Decomposer::updateContextVector(JADE *optimizer)
{
    vector<double> v = optimizer->getCollaborator();
    double newBestCandidate = optimizer->calculateFitnessValue(v);

    if (newBestCandidate BETTER_THAN bestAchievedFitness)
    {
        for (unsigned ld = 0; ld<v.size(); ld++)
            contextVector[optimizer->coordinates[ld]] = v[ld];

        bestAchievedFitness = newBestCandidate;
    }
}



void Decomposer::buildContextVector()
{
    for (unsigned j = 0; j<optimizers.size(); ++j)
    {
        vector<double> v = optimizers[j]->getCollaborator();
        double newBestCandidate = optimizers[j]->calculateFitnessValue(v);

        if (newBestCandidate BETTER_THAN bestAchievedFitness)
        {
            for (unsigned ld = 0; ld<v.size(); ld++)
                contextVector[optimizers[j]->coordinates[ld]] = v[ld];
            bestAchievedFitness = newBestCandidate;
        }
    }
    //bestAchievedFitness = costFunction(contextVector);
}



void Decomposer::randomGrouping()
{
    if ( optimizers.size() && this->applyRandomGrouping )
    {
        shuffle(coordinates.begin(), coordinates.end(), eng);
        //setOptimizersCoordinatesAndEvaluatePopulation();

        unsigned numOfCoordinatesPerSubgroup = coordinates.size() / optimizers.size();
        for (unsigned i = 0; i < optimizers.size(); ++i)
        {
            optimizers[i]->setCoordinates(&(coordinates[i*numOfCoordinatesPerSubgroup]), numOfCoordinatesPerSubgroup);
            for (int k = 0; k < sizeOfSubcomponents; ++k)
                coordinateToOptimizer[coordinates[i*numOfCoordinatesPerSubgroup + k]] = optimizers[i];
            optimizers[i]->loadIndividuals(population);
            optimizers[i]->evaluateParents();
            optimizers[i]->updateIndexOfBest();
        }

    }
}


void Decomposer::sortPopulation(typeOfSort t)
{
    popSort.clear();
    fitnessValues.clear();
    for (unsigned i = 0; i < population.size(); ++i)
    {
        fitnessValues.push_back(CCOptimizer.computeFitness(population[i]));
        popSort.push_back(i);
    }
    sort(popSort.begin(), popSort.end(), Decomposer::doCompareIndividuals(fitnessValues, t));
}



