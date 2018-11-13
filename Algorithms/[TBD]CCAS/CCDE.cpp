//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : CCDE.cpp
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

#include "CCDE.h"
#include "Decomposer.h"
#include "JADE.h"



//******************************************************************************************/
//
//
//
//******************************************************************************************/
CCDE::CCDE()
{
    JADE_c = 0.1;
    JADE_p = 0.1;

    JADE_mutationStrategy = 2;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer * CCDE::createDecomposer(unsigned sizeOfSubcomponents, unsigned individualsPerSubcomponent, bool RG)
{
    vector<unsigned> allCoordinates;
    for(unsigned i=0; i<problemDimension; ++i)
        allCoordinates.push_back(i);

    if ( RG )
        shuffle(allCoordinates.begin(), allCoordinates.end(), eng);

    unsigned seed = unifRandom(eng)*100000;

    return new Decomposer(*this, seed, allCoordinates, sizeOfSubcomponents, individualsPerSubcomponent, population, contextVector, RG);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer * CCDE::createDecomposer(unsigned indexOfDecomposerCharacteristics, bool RG)
{
    vector<unsigned> allCoordinates;
    for (unsigned i = 0; i<problemDimension; ++i)
        allCoordinates.push_back(i);

    if (RG)
        shuffle(allCoordinates.begin(), allCoordinates.end(), eng);

    unsigned sizeOfSubcomponents = subcomponentSizes[indexOfDecomposerCharacteristics];
    unsigned individualsPerSubcomponent = numIndividualsPerSubcomponents[indexOfDecomposerCharacteristics];

    unsigned seed = unifRandom(eng) * 100000;

    Decomposer *dec = new Decomposer(*this, seed, allCoordinates, sizeOfSubcomponents, individualsPerSubcomponent, population, contextVector, RG);
    dec->indexOfDecomposerCharacteristics = indexOfDecomposerCharacteristics;
    return dec;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::createSetOfDecomposers(vector<unsigned> &sizesOfSubcomponents, vector<unsigned> &individualsPerSubcomponent, bool RG)
{
    vector<unsigned> allCoordinates;
    for (unsigned i = 0; i<problemDimension; ++i)
        allCoordinates.push_back(i);

    if ( RG )
        shuffle(allCoordinates.begin(), allCoordinates.end(), eng);

    for (unsigned s = 0; s < decomposers.size(); ++s)
        delete decomposers[s];
    decomposers.clear();

    unsigned seed = unifRandom(eng) * 100000;

    for (unsigned k = 0; k < sizesOfSubcomponents.size(); ++k)
        decomposers.push_back(new Decomposer(*this, seed, allCoordinates, sizesOfSubcomponents[k], individualsPerSubcomponent[k], population, contextVector, RG));
}





//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::initPopulation(unsigned numOfIndividuals)
{
    for (unsigned i = 0; i<numOfIndividuals; ++i)
    {
        vector< double > position;
        for (unsigned int d = 0; d<problemDimension; ++d)
            position.push_back(lowerLimit + unifRandom(eng)*(upperLimit - lowerLimit));
        population.push_back(position);
        fitnessValues.push_back(0);
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
double CCDE::computeFitness(vector<double> &x)
{
    numberOfEvaluations++;
    return fitness->compute(&x[0]);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
CCDE::~CCDE()
{
    for (unsigned i = 0; i < decomposers.size(); ++i)
        delete decomposers[i];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::initContextVector()
{
    contextVector.resize(problemDimension);

    ///Initialize the context vector (just pick the first individual, which is randomly generated)
    for (unsigned i = 0; i < problemDimension; ++i)
        contextVector[i] = population[0][i];

    globalBestFitness = computeFitness(contextVector);

    numberOfEvaluations++;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
double CCDE::getFinalFitnessValue()
{
    return globalBestFitness;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimize(Fitness* _function, unsigned int maxNumberOfEvaluations,
                    unsigned _sizeOfSubcomponents,
                    unsigned individualsPerSubcomponent,
                    unsigned seed, unsigned nGenPerCycle)

{
    cout << "Using standard CC optimizer..." << endl;

    eng.seed(seed);
    numberOfEvaluations = 0;

    fitness = _function;

    clock_t startTime = clock();

    problemDimension = fitness->getDimension();

    lowerLimit = fitness->getMinX();

    upperLimit = fitness->getMaxX();

    initPopulation(individualsPerSubcomponent);

    initContextVector();

    Decomposer *dec = createDecomposer(_sizeOfSubcomponents, individualsPerSubcomponent, true);

    ostringstream funcId;
    funcId << _function->getID();
    string outFileName = "STDCC_Function_f" + funcId.str() + ".csv";
    ofstream outFile;
    outFile.open(outFileName);
    outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << _sizeOfSubcomponents << ";" << individualsPerSubcomponent << ";\n";

    for (unsigned ite = 0; numberOfEvaluations<maxNumberOfEvaluations; ++ite)
    {
        for(unsigned  j = 0; j < dec->optimizers.size(); ++j)
        {
            JADE *optimizer = dec->optimizers[j];
            optimizer->loadIndividuals(dec->population);

            optimizer->optimize(nGenPerCycle);

            optimizer->storeIndividuals(dec->population);
        }

        dec->buildContextVector();
        if ( dec->bestAchievedFitness < globalBestFitness )
            globalBestFitness = dec->bestAchievedFitness;

        dec->randomGrouping();

        outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << _sizeOfSubcomponents << ";" << individualsPerSubcomponent << ";\n";
        cout << numberOfEvaluations << "   " << globalBestFitness << "    " << _sizeOfSubcomponents << "    " << individualsPerSubcomponent << "\n";
    }
    clock_t stopTime = clock();
    elapsedTime = ((double)(stopTime - startTime))/CLOCKS_PER_SEC;
    outFile.close();
    delete dec;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
unsigned CCDE::findDecomposerWithMaxValueFunction()
{
    double bestValue = -1.0E20;
    int bestDecomposer = 0;
    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        if (decomposers[s]->valueFunction > bestValue && fabs(decomposers[s]->valueFunction)>0)
        {
            bestValue = decomposers[s]->valueFunction;
            bestDecomposer = s;
        }
    }
    return bestDecomposer;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::copySearchState(unsigned bestFitnessDecomposer, unsigned activeDecomposer)
{
    decomposers[bestFitnessDecomposer]->sortPopulation(Decomposer::fromBestToWorst);

    if (decomposers[bestFitnessDecomposer]->population.size() < decomposers[activeDecomposer]->population.size())
    {
        decomposers[activeDecomposer]->sortPopulation(Decomposer::fromWorstToBest);

        for (unsigned i = 0; i < decomposers[bestFitnessDecomposer]->population.size(); ++i)
        {
            unsigned ii = decomposers[bestFitnessDecomposer]->popSort[i];
            unsigned jj = decomposers[activeDecomposer]->popSort[i];
            for (unsigned d = 0; d < problemDimension; ++d)
                decomposers[activeDecomposer]->population[jj][d] = decomposers[bestFitnessDecomposer]->population[ii][d];
        }
    }
    else
    {
        for (unsigned i = 0; i < decomposers[activeDecomposer]->population.size(); ++i)
        {
            unsigned ii = decomposers[bestFitnessDecomposer]->popSort[i];
            for (unsigned d = 0; d < problemDimension; ++d)
                decomposers[activeDecomposer]->population[i][d] = decomposers[bestFitnessDecomposer]->population[ii][d];
        }
    }

    int ii = decomposers[activeDecomposer]->population.size() - 1;
    for (unsigned d = 0; d < problemDimension; ++d)
        decomposers[activeDecomposer]->population[ii][d] = decomposers[bestFitnessDecomposer]->contextVector[d];

    decomposers[activeDecomposer]->contextVector = decomposers[bestFitnessDecomposer]->contextVector;

    for (unsigned j = 0; j < decomposers[activeDecomposer]->optimizers.size(); ++j)
    {
        JADE *optimizer = decomposers[activeDecomposer]->optimizers[j];
        optimizer->loadIndividuals(decomposers[activeDecomposer]->population);
        numberOfEvaluations += optimizer->evaluateParents();
        optimizer->updateIndexOfBest();
    }

}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::selectNextPopulation(unsigned activeDecomposer, unsigned bestFitnessDecomposer)
{
    vector< vector<double> > population;

    decomposers[activeDecomposer]->sortPopulation(Decomposer::fromBestToWorst);
    decomposers[bestFitnessDecomposer]->sortPopulation(Decomposer::fromBestToWorst);

    for (unsigned i = 0; i < decomposers.size(); ++i)
        population.push_back(decomposers[i]->contextVector);

    int i = 0;
    while (population.size() < decomposers[activeDecomposer]->population.size() &&
            i < decomposers[bestFitnessDecomposer]->population.size())
        population.push_back(decomposers[bestFitnessDecomposer]->population[decomposers[bestFitnessDecomposer]->popSort[i++]]);

    while (population.size() < decomposers[activeDecomposer]->population.size())
        population.push_back(decomposers[activeDecomposer]->population[decomposers[activeDecomposer]->popSort[i++]]);

    decomposers[activeDecomposer]->population = population;
    decomposers[activeDecomposer]->contextVector = decomposers[bestFitnessDecomposer]->contextVector;
    decomposers[activeDecomposer]->setOptimizersCoordinatesAndEvaluatePopulation();
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::selectNextPopulation(unsigned activeDecomposer)
{
    vector< vector<double> > population;

    decomposers[activeDecomposer]->sortPopulation(Decomposer::fromBestToWorst);
    for (unsigned i = 0; i < decomposers.size(); ++i)
        population.push_back(decomposers[i]->contextVector);

    int i = 0;
    while ( population.size() < decomposers[activeDecomposer]->population.size() )
        population.push_back(decomposers[activeDecomposer]->population[decomposers[activeDecomposer]->popSort[i++]]);

    decomposers[activeDecomposer]->population = population;
    decomposers[activeDecomposer]->setOptimizersCoordinatesAndEvaluatePopulation();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/

void CCDE::broadcastSearchState(unsigned sourceDecomposer, double randomRatio)
{
    decomposers[sourceDecomposer]->sortPopulation(Decomposer::fromBestToWorst);

    //homogenize JADE parameters among optimizers and decomposers (not really necessary)
    double JADE_mu_cr = 0, JADE_mu_ff = 0;
    for (int i = 0; i < decomposers[sourceDecomposer]->optimizers.size(); ++i)
    {
        JADE_mu_cr += decomposers[sourceDecomposer]->optimizers[i]->JADE_mu_cr;
        JADE_mu_ff += decomposers[sourceDecomposer]->optimizers[i]->JADE_mu_ff;
    }
    JADE_mu_cr /= decomposers[sourceDecomposer]->optimizers.size();
    JADE_mu_ff /= decomposers[sourceDecomposer]->optimizers.size();

    for (int i = 0; i < decomposers[sourceDecomposer]->optimizers.size(); ++i)
    {
        decomposers[sourceDecomposer]->optimizers[i]->JADE_mu_cr = JADE_mu_cr;
        decomposers[sourceDecomposer]->optimizers[i]->JADE_mu_ff = JADE_mu_ff;
    }


    //let all decomposers start with the same random seed
    unsigned seed = unifRandom(eng) * 100000;
    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        decomposers[s]->eng.seed(seed);
    }

    //create a buffer of random individuals
    vector< vector<double> > sourcePopulation;
    for (int i = 0; i < maxPopSize; ++i)
    {
        vector<double> v;
        for (unsigned d = 0; d < problemDimension; ++d)
            v.push_back( lowerLimit + unifRandom(eng) * (upperLimit - lowerLimit) );
        sourcePopulation.push_back(v);
    }


    //copy the source population
    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        if (s != sourceDecomposer)
        {
            decomposers[s]->population[0] = decomposers[sourceDecomposer]->contextVector;

            for (unsigned i = 1; i < decomposers[s]->population.size(); ++i)
                if ( i < decomposers[sourceDecomposer]->population.size() )
                {
                    int ii = decomposers[sourceDecomposer]->popSort[i];

                    for (unsigned d = 0; d < problemDimension; ++d)
                        decomposers[s]->population[i][d] = decomposers[sourceDecomposer]->population[ii][d];
                }

            decomposers[s]->contextVector = decomposers[sourceDecomposer]->contextVector;
            decomposers[s]->bestAchievedFitness = decomposers[sourceDecomposer]->bestAchievedFitness;
            decomposers[s]->prevBestFitness = decomposers[sourceDecomposer]->bestAchievedFitness;
            decomposers[s]->coordinates = decomposers[sourceDecomposer]->coordinates;
        }

        for (int i = 0; i < decomposers[s]->optimizers.size(); ++i)
        {
            decomposers[s]->optimizers[i]->JADE_mu_cr = JADE_mu_cr;
            decomposers[s]->optimizers[i]->JADE_mu_ff = JADE_mu_ff;
        }
    }

    //overwrite with some random individuals
    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        for (unsigned i = decomposers[s]->population.size()*randomRatio; i < decomposers[s]->population.size(); ++i)
        {
            decomposers[s]->population[i] = sourcePopulation[i];
        }
    }

    sourcePopulation.clear();
    for (int i = 0; i < decomposers[sourceDecomposer]->population.size(); ++i)
        sourcePopulation.push_back(decomposers[sourceDecomposer]->population[decomposers[sourceDecomposer]->popSort[i]]);

    decomposers[sourceDecomposer]->population = sourcePopulation;

    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        decomposers[s]->setOptimizersCoordinatesAndEvaluatePopulation();
    }

}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimizeSubcomponents(Decomposer *dec, unsigned nGenPerIteration)
{
    dec->functionEvaluations = 0;

    for (unsigned j = 0; j < dec->optimizers.size(); ++j)
    {
        JADE *optimizer = dec->optimizers[j];
        optimizer->loadIndividuals(dec->population);
        dec->functionEvaluations += optimizer->optimize(nGenPerIteration);
        optimizer->storeIndividuals(dec->population);
    }
    dec->buildContextVector();
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimize_CCAS(Fitness* _function,
                         unsigned _maxNumberOfEvaluations,
                         vector<unsigned> &subcomponentSizes,
                         vector<unsigned> &numIndividualsPerSubcomponents,
                         unsigned seed,
                         unsigned nGenPerCycle,
                         unsigned LCP,
                         double rho,
                         double lambda,
                         double delta,
                         unsigned sigma)
{
    cout << "Using CCAS..." << endl;

    maxNumberOfEvaluations = _maxNumberOfEvaluations;
    unsigned fitnessEvaluationsBetweenLearning = maxNumberOfEvaluations*delta;
    unsigned maxFitnessEvaluationsForLearning = maxNumberOfEvaluations*lambda;
    maxPopSize = *max_element(numIndividualsPerSubcomponents.begin(), numIndividualsPerSubcomponents.end());
    eng.seed(seed);
    numberOfEvaluations = 0;
    fitness = _function;
    clock_t startTime = clock();
    problemDimension = fitness->getDimension();
    lowerLimit = fitness->getMinX();
    upperLimit = fitness->getMaxX();

    initPopulation(maxPopSize);

    initContextVector();

    createSetOfDecomposers(subcomponentSizes, numIndividualsPerSubcomponents, true);

    unsigned stagnationFlag = 0;
    unsigned totalLearningCost = 0;
    unsigned activeDecomposer = 0;
    unsigned numberOfEvaluationsFromLastLearning = 0;

    cout << "Number of generations per cycle =" << nGenPerCycle << endl;
    cout << "LCP =" << LCP << endl;
    cout << "rho =" << rho << endl;
    cout << "Max fitness evaluations for comparison phases =" << maxFitnessEvaluationsForLearning << endl;
    cout << "Learning every =" << fitnessEvaluationsBetweenLearning << " function evaluations" << endl;
    cout << "stagnation limit =" << sigma << endl << endl;

    ostringstream funcId;
    funcId << _function->getID();
    string outFileName = "CCAS_Function_f" + funcId.str() + ".csv";
    ofstream outFile;
    outFile.open(outFileName);
    outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << ";" << ";\n";

    for (unsigned s = 0; s < decomposers.size(); ++s)
    {
        decomposers[s]->bestAchievedFitness = globalBestFitness;
    }

    double prevBestFitness = globalBestFitness;
    bool comparisonPhase = true;
    for (unsigned ite = 0; numberOfEvaluations<maxNumberOfEvaluations; ++ite)
    {
        if (comparisonPhase)
        {
            double ff = globalBestFitness;
            for (unsigned s = 0; s < decomposers.size(); ++s)
            {
                decomposers[s]->learningCost = 0;

                for (unsigned learningCycle = 0; learningCycle<LCP; ++learningCycle)
                {
                    //Iterate the search algorithm in each subcomponent
                    optimizeSubcomponents(decomposers[s], nGenPerCycle);

                    decomposers[s]->learningCost += decomposers[s]->functionEvaluations;
                    totalLearningCost += decomposers[s]->functionEvaluations;

                    decomposers[s]->randomGrouping();

                    ff = min(ff, decomposers[s]->bestAchievedFitness);
                    outFile << numberOfEvaluations << ";" << ff << ";" << decomposers[activeDecomposer]->sizeOfSubcomponents << ";" << decomposers[activeDecomposer]->individualsPerSubcomponent << ";\n";
                }

                decomposers[s]->valueFunction = max(0.0, globalBestFitness - decomposers[s]->bestAchievedFitness) / decomposers[s]->learningCost;
            }


            comparisonPhase = false;

            sort(decomposers.begin(), decomposers.end(), doCompareDecomposers());

            cout << "Ranking of configurations:" << endl;
            vector<Decomposer*>::iterator decIte = decomposers.begin();
            for (; decIte != decomposers.end(); ++decIte)
                cout << (*decIte)->sizeOfSubcomponents << "/" << (*decIte)->individualsPerSubcomponent << " -> " << (*decIte)->valueFunction << endl;

            activeDecomposer = 0;

            //Find the decomposer with the best fitness
            unsigned bestFitnessDecomposer = 0;
            double bestFitness = decomposers[0]->bestAchievedFitness;
            for (unsigned s = 1; s < decomposers.size(); ++s)
                if (decomposers[s]->bestAchievedFitness < bestFitness)
                {
                    bestFitness = decomposers[s]->bestAchievedFitness;
                    bestFitnessDecomposer = s;
                }
            selectNextPopulation(activeDecomposer, bestFitnessDecomposer);

            cout << "Selected configuration: subcomponents of size " << decomposers[activeDecomposer]->sizeOfSubcomponents << " with " << decomposers[activeDecomposer]->individualsPerSubcomponent << " individuals" << endl;
            cout << ite << " " << "NOE=" << numberOfEvaluations << " " << std::scientific << bestFitness << endl;
            cout << "== End of comparison phase ==" << endl;

            numberOfEvaluationsFromLastLearning = 0;

            globalBestFitness = bestFitness;

            ///
            decomposers[activeDecomposer]->prevBestFitness = decomposers[activeDecomposer]->bestAchievedFitness;
            ///
        }
        else
        {
            //Iterate the search algorithm in each subcomponent of the selected decomposer
            optimizeSubcomponents(decomposers[activeDecomposer], nGenPerCycle);

            numberOfEvaluationsFromLastLearning += decomposers[activeDecomposer]->functionEvaluations;

            //Update best fitness of the s-th decomposer
            globalBestFitness = decomposers[activeDecomposer]->bestAchievedFitness;

            decomposers[activeDecomposer]->randomGrouping();

            cout << ite << " " << "NOE=" << numberOfEvaluations << " " << std::scientific << globalBestFitness << endl;
            if (fabs(decomposers[activeDecomposer]->prevBestFitness) > 0)
            {
                if (fabs(decomposers[activeDecomposer]->bestAchievedFitness - decomposers[activeDecomposer]->prevBestFitness) / decomposers[activeDecomposer]->prevBestFitness < 1.0E-06)
                    stagnationFlag++;
                else
                    stagnationFlag = 0;
            }
            //Check if learning should be reactivated
            if (((totalLearningCost<maxFitnessEvaluationsForLearning && numberOfEvaluationsFromLastLearning >= fitnessEvaluationsBetweenLearning) || stagnationFlag>sigma) && decomposers.size()>1)
            {
                cout << "=== Reactivating comparison phase ===" << endl;
                comparisonPhase = true;
                stagnationFlag = 0;
                broadcastSearchState(activeDecomposer, 1-rho);
                prevBestFitness = globalBestFitness;
            }

            decomposers[activeDecomposer]->prevBestFitness = decomposers[activeDecomposer]->bestAchievedFitness;
        }

        outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << decomposers[activeDecomposer]->sizeOfSubcomponents << ";" << decomposers[activeDecomposer]->individualsPerSubcomponent << ";\n";
    }
    clock_t stopTime = clock();
    outFile.close();
    elapsedTime = ((double)(stopTime - startTime)) / CLOCKS_PER_SEC;
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
unsigned CCDE::selectDecomposerSoftMax(double tau, vector<double> &valueFunction)
{
    unsigned numOfDecomposers = valueFunction.size();
    vector<double> decomposersProb(numOfDecomposers, 0);

    //double norm = sqrt(inner_product(valueFunction.begin(), valueFunction.end(), valueFunction.begin(), 0.0));
    double norm = *max_element(valueFunction.begin(), valueFunction.end());
    if (norm < 0.00001)
        return unifRandom(eng)*numOfDecomposers;
    transform(valueFunction.begin(), valueFunction.end(), valueFunction.begin(), [norm](double d) -> double { return d / norm; });

    double den = 0.0;
    for (unsigned i = 0; i < numOfDecomposers; ++i)
    {
        decomposersProb[i] = exp(valueFunction[i] / tau);
        den += decomposersProb[i];
    }
    transform(decomposersProb.begin(), decomposersProb.end(), decomposersProb.begin(), [den](double d) -> double { return d / den; });

    double p = unifRandom(eng);
    double cumProb = 0;
    for (unsigned i = 0; i < numOfDecomposers; ++i)
    {
        cumProb += decomposersProb[i];
        if (cumProb >= p)
            return i;
    }
    return numOfDecomposers - 1;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
unsigned CCDE::selectDecomposerSoftMax(double tau, vector<Decomposer*> &dec)
{
    unsigned numOfDecomposers = dec.size();

    if (numOfDecomposers == 0)
    {
        cerr << "No decomposers available" << endl;
        exit(1);
    }

    vector<double> valueFunction(numOfDecomposers, 0);

    for (unsigned i = 0; i < numOfDecomposers; ++i)
        valueFunction[i] = dec[i]->valueFunction;

    return selectDecomposerSoftMax(tau, valueFunction);
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimize_MLSOFT(Fitness* _function,
                           unsigned int maxNumberOfEvaluations,
                           vector<unsigned> &_subcomponentSizes,
                           unsigned numIndividualsPerSubcomponents,
                           unsigned seed,
                           double tau,
                           unsigned nGenPerCycle)
{
    cout << "Using MLSoft..." << endl;
    cout << "tau = " << tau << endl;
    cout << "Individuals per subcomponent = " << numIndividualsPerSubcomponents << endl<<endl;

    subcomponentSizes = _subcomponentSizes;
    eng.seed(seed);
    numberOfEvaluations = 0;
    unsigned currentDecomposer = unifRandom(eng)*subcomponentSizes.size();
    unsigned prevDecomposer = currentDecomposer;
    unsigned numOfDecomposers = subcomponentSizes.size();
    vector<double> valueFunction(numOfDecomposers, 0);
    vector<unsigned> decomposersCounters(numOfDecomposers, 0);

    fitness = _function;

    clock_t startTime = clock();

    problemDimension = fitness->getDimension();

    lowerLimit = fitness->getMinX();

    upperLimit = fitness->getMaxX();

    initPopulation(numIndividualsPerSubcomponents);

    initContextVector();

    Decomposer *dec = createDecomposer(subcomponentSizes[currentDecomposer], numIndividualsPerSubcomponents, true);

    ostringstream funcId;
    funcId << _function->getID();
    string outFileName = "MLSoft_Function_f" + funcId.str() + ".csv";
    ofstream outFile;
    outFile.open(outFileName);
    outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << subcomponentSizes[currentDecomposer] << ";" << numIndividualsPerSubcomponents << ";\n";

    double prevBestFitness = 0;

    //main cycle
    for (unsigned ite = 0; numberOfEvaluations<maxNumberOfEvaluations; ++ite)
    {

        dec->randomGrouping();

        //choose the decomposer
        currentDecomposer = selectDecomposerSoftMax(tau, valueFunction);

        decomposersCounters[currentDecomposer]++;

        //set the chosen decomposer
        if (currentDecomposer != prevDecomposer)
        {
            for (unsigned i = 0; i < population.size(); ++i)
                fitnessValues[i] = computeFitness(population[i]);

            dec->setSubcomponentsOfEqualSize(subcomponentSizes[currentDecomposer], fitnessValues);

            prevBestFitness = globalBestFitness;
            prevDecomposer = currentDecomposer;
        }


        //Iterate the search algorithm in each component
        unsigned lastNumberOfEvaluations = 0;
        for(unsigned  j = 0; j < dec->optimizers.size(); ++j)
        {
            //cout << "optimize subcomponent" << endl;
            JADE *optimizer = dec->optimizers[j];

            optimizer->loadIndividuals(population);
            optimizer->optimize(nGenPerCycle);
            optimizer->storeIndividuals(population);
        }

        dec->buildContextVector();
        if ( dec->bestAchievedFitness < globalBestFitness)
            globalBestFitness = dec->bestAchievedFitness;


        double r = 0;
        if ( ite > 0 )
        {
            //reward
            r = (prevBestFitness - globalBestFitness) / fabs(prevBestFitness);

            //update value function
            valueFunction[currentDecomposer] = (valueFunction[currentDecomposer] * decomposersCounters[currentDecomposer] + r) / (decomposersCounters[currentDecomposer] + 1);
        }

        prevBestFitness = globalBestFitness;

        outFile << numberOfEvaluations << ";" << globalBestFitness << ";" << subcomponentSizes[currentDecomposer] << ";" << numIndividualsPerSubcomponents << ";\n";

        cout << numberOfEvaluations << "   " << globalBestFitness << "   " << subcomponentSizes[currentDecomposer] << "   " << numIndividualsPerSubcomponents << "\n";
    }
    clock_t stopTime = clock();
    elapsedTime = ((double)(stopTime - startTime)) / CLOCKS_PER_SEC;
    delete dec;
    outFile.close();
}
