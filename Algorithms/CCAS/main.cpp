//=======================================================================================
// Cooperative Coevolution with Adaptive Subcomponents
//=======================================================================================
// Name        : main.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Pawel Topa
//               Jaroslaw Was
// Version     : v1.0
// Created on  : Mar 20, 2016
//
// More details on the following paper:
//
// Trunfio, G.A., Topa, P., Was, J. 'A New Algorithm for Adapting the Configuration
// of Subcomponents in Large-Scale Optimization with Cooperative Coevolution, submitted'
//
//=======================================================================================


#define _USE_MATH_DEFINES
#include "CCDE.h"
#include "Benchmarks2010.h"
#include "Header.h"


vector<unsigned> subcomponentSizes = { 2, 5, 10, 25, 50, 100, 200, 500, 1000 };  //hard coded
vector<unsigned> numIndividualsPerSubcomponents = { 20, 20, 20, 30, 50, 50, 100, 100, 200 };  //hard coded
unsigned seed = 1;
unsigned type = 2;
unsigned int numberOfEvaluations = 3000000;
unsigned functionIndex = 1;
double tau = 0.5;
unsigned sizeOfSubcomponents = 100;
unsigned numOfIndividuals = 50;
unsigned nGenPerCycle = 4;
unsigned LCP = 2;
double rho = 0.2;
double lambda = 0.15;
double delta = 0.05;
double sigma = 10;


// create new object of class with default setting
Fitness* generateFuncObj(int funcID)
{
    Fitness *fp=NULL;
    using namespace CEC2010;
    if (funcID == 1)       fp = new F1();
    else if (funcID == 2)  fp = new F2();
    else if (funcID == 3)  fp = new F3();
    else if (funcID == 4)  fp = new F4();
    else if (funcID == 5)  fp = new F5();
    else if (funcID == 6)  fp = new F6();
    else if (funcID == 7)  fp = new F7();
    else if (funcID == 8)  fp = new F8();
    else if (funcID == 9)  fp = new F9();
    else if (funcID == 10) fp = new F10();
    else if (funcID == 11) fp = new F11();
    else if (funcID == 12) fp = new F12();
    else if (funcID == 13) fp = new F13();
    else if (funcID == 14) fp = new F14();
    else if (funcID == 15) fp = new F15();
    else if (funcID == 16) fp = new F16();
    else if (funcID == 17) fp = new F17();
    else if (funcID == 18) fp = new F18();
    else if (funcID == 19) fp = new F19();
    else if (funcID == 20) fp = new F20();
    else
    {
        cerr << "Fail to locate Specified Function Index" << endl;
        exit(-1);
    }

    return fp;
}


//================================================================================================
//
//
//================================================================================================
void help()
{
    cout << "Cooperative Coevolution with Adaptive Subcomponents" << endl;
    cout << "More details on the following paper:" << endl;
    cout << "Trunfio, G.A., Topa, P., Was, J. 'An Algorithm for Adapting the Configuration"  << endl;
    cout << "of Subcomponents in Large-Scale Optimization with Cooperative Coevolution', submitted" << endl;
    cout << "trunfio@uniss.it" << endl << endl;


    cout << "Syntax:" << endl;
    cout << "CCAS -f x  [-h|-?] [-NFE x][-type x][-seed x][-tau x][-sizeOfSubs x][-numOfInd x][-nGenPC x][-LCP x][-rho x][-lambda x][-delta x][-sigma x]" << endl;
    cout << "Parameters:  " << endl;
    cout << "  -f x             -> optimize function x from CEC 2010 LSGO test bed with in [1, 20]" << endl;
    cout << "  -h|-?            -> command syntax" << endl;
    cout << "  -NFE x           -> set x as the maximum number of fitness evaluations (default is " << numberOfEvaluations << ")"  << endl;
    cout << "  -type t          -> set t as type of optimization: t=0 for standard CC, t=1 for MLSoft, t=2 for CCAS (default is " << type << ")" << endl;
    cout << "  -seed x          -> set x as random seed (default is time(0))" << endl;
    cout << "  -tau x           -> set the tau parameter for MLSoft (delault is " << endl;
    cout << "  -sizeOfSubs x    -> set the size of subcomponents for the standard CC" << endl;
    cout << "  -numOfInd x      -> set the number of individuals for standard CC/MLSoft" << endl;
    cout << "  -nGenPC x        -> set the number of generations per cycle" << endl;
    cout << "  -LCP x           -> set the length in cycles for the comparison phases of CCAS" << endl;
    cout << "  -rho x           -> set the rho parameter for CCAS (see the article) (default is " << rho << ")" << endl;
    cout << "  -lambda x        -> set the lambda parameter for CCAS (see the article) (default is " << lambda << ")" << endl;
    cout << "  -delta x         -> set the delta parameter for CCAS (see the article) (default is " << delta << ")" << endl;
    cout << "  -sigma x         -> set the sigma parameter for CCAS (see the article) (default is " << sigma << ")" << endl;
}


//================================================================================================
//
//
//================================================================================================
void parseParameters(int argc, char* argv[])
{
    for (int i = 1; i<argc; i++)
    {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-?") == 0))
        {
            help();
            exit(0);
        }
        else if (strcmp(argv[i], "-f") == 0)         // test function index
        {
            functionIndex = atoi(argv[++i]);
            if (functionIndex < 1 || functionIndex>20)
            {
                cerr << "Invalid function index" << endl;
                help();
                exit(0);
            }
        }
        else if (strcmp(argv[i], "-NFE") == 0)       // max number of function evaluations
        {
            numberOfEvaluations = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-type") == 0)
        {
            type = atoi(argv[++i]);
            if (type < 0 || type>2)
            {
                cerr << "Invalid type of algorithm" << endl;
                help();
                exit(0);
            }
        }
        else if (strcmp(argv[i], "-seed") == 0)
        {
            seed = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-tau") == 0)
        {
            tau = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-sizeOfSubs") == 0)
        {
            sizeOfSubcomponents = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-numOfInd") == 0)
        {
            numOfIndividuals = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-nGenPC") == 0)
        {
            nGenPerCycle = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-LCP") == 0)
        {
            LCP = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-rho") == 0)
        {
            rho = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-lambda") == 0)
        {
            lambda = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-delta") == 0)
        {
            delta = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-sigma") == 0)
        {
            sigma = atoi(argv[++i]);
        }
        else
        {
            cerr << "Fatal error: invalid parameter: " << argv[i] << endl;
            help();
            exit(1);
        }
    }
}


void main(int argc, char* argv[])
{
    seed = time(0);

    parseParameters(argc, argv);

    if (functionIndex == -1)
    {
        help();
        exit(1);
    }

    Fitness* f = generateFuncObj(functionIndex);

    if (f == NULL)
    {
        cerr << "Unable to create fitness function" << endl;
        exit(1);
    }

    cout << "Optimizing f" << functionIndex << " from CEC 2010 LSGO test suite" << endl;

    CCDE ccde;

    if ( type==2 )
        ccde.optimize_CCAS(f, numberOfEvaluations, subcomponentSizes, numIndividualsPerSubcomponents, seed, nGenPerCycle, LCP, rho, lambda, delta, sigma);
    else if (type == 1)
        ccde.optimize_MLSOFT(f, numberOfEvaluations, subcomponentSizes, numOfIndividuals, seed, tau, nGenPerCycle);
    else if (type == 0)
        ccde.optimize(f, numberOfEvaluations, sizeOfSubcomponents, numOfIndividuals, seed, nGenPerCycle);

}
