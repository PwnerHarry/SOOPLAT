/*
 * =====================================================================================
 *
 *       Filename:  CCVIL.h
 *
 *
 *        Version:  1.0
 *        Created:  02/24/2011 08:37:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wenxiang Chen (http://cs-chen.net), chenwx.ustc@gmail.com
 *        Company:  Nature Inspired Computation and Application Laboratory (NICAL), USTC
 *
 * =====================================================================================
 */

#ifndef _CCVIL_H
#define _CCVIL_H

#include	<EALib/PopulationT.h>
#include	<EALib/IndividualT.h>
#include	<Rng/Normal.h>
#include	<Rng/Uniform.h>
#include	<vector>
#include	<iostream>
#include	<algorithm>
#include	<climits>
#include	<sys/stat.h>
#include	<sys/types.h>
#include	<sys/time.h>
#include	<fcntl.h>
#include	<string>
#include	<sstream>
#include	<ctime>
#include	<cfloat>
#include	<algorithm>


#include "benchmark/Benchmarks.h"
#include "Archive.h"

using namespace std;

class CCVIL{
	protected:
		vector< PopulationT<double> > pop;
		vector< vector<unsigned> > groupInfo;
		vector<unsigned> samplingPoints;

		IndividualT<double>* bestCand;
		Benchmarks* fp;
		unsigned* p; // random permutation
		unsigned* lookUpGroup;
		unsigned MaxFitEval; // total fitness evaluation limitation for optimization
		unsigned curFitEval; // current fitness evaluation limitation for optimization
		unsigned cycle; 

		unsigned* randPerm(unsigned N);
		unsigned lowerThreshold;
		unsigned upperThreshold;
		RunParameter* param;
		double lastBestValue;
		double* groupCR;
		double* groupF;
		unsigned* failCounter;
		unsigned long fes;
		double bestFit;
		bool* impreciseGroup;
		unsigned bestCandIndex;

		// in result folder
		FILE *resultFP;
		vector<double> resultRec;
		FILE *timeFP;
		vector<double> timeRec;

		// in trace folder
		FILE *groupFP; 
		vector<unsigned> groupRec;
		FILE *groupFesFP; 
		vector<unsigned> groupFesRec;
		FILE *fesFP;
		vector<unsigned> fesRec;
		FILE *valFP;
		vector<double> valRec;

		void learningStage();
		void optimizationStage();
		bool sameGroup(unsigned v1, unsigned v2);
		unsigned JADECC(unsigned index, bool learnStageFlag);
		void captureInter(unsigned curDim, unsigned lastDim);
		void popInit();
		void popInitZeros ( );
		void popGenerate(bool learnStageFlag);
		void popSizeVary ( double factor );
		void randFCR(unsigned NP, double CRm, double CRsigma, double Fm, double Fsigma, double *&F, double *&CR);
		double cauchyRnd(double mu, double delta);
		void printArray(double* a, unsigned D);
		void printArray(unsigned* a, unsigned D);
		void printPopulation(PopulationT<double> printPop);
		void printPopulation(IndividualT<double> printIndiv);
		void print2Dvector(vector< PopulationT<double> > vector2D);
		void print2Dvector(vector< vector<unsigned> > vector2D);
		void printFitness(PopulationT<double> printPop);
		void printVector(vector<unsigned> v);
		void printVector(vector<double> v);
		void printVector(vector<bool> v);
		double sum(vector<double> vec);
		unsigned sum ( unsigned* arr, unsigned N );
		double mean ( double* arr, unsigned size );
		vector<double> dotMultiply(vector<double> v1, vector<double> v2);
		void gnR1R2(unsigned NP1, unsigned NP2, unsigned *r1, unsigned *r2);
		PopulationT<double> combinePopulation(PopulationT<double> p1, PopulationT<double> p2);
		void findPbestIndex(PopulationT<double> inPop, unsigned pNP, unsigned* indBest);
		void boundConstrain(PopulationT<double> &vi, PopulationT<double> offsprings, int LB, int UB, vector<unsigned> vecIndex);
		string itos ( int i );
		void sampleInfo ( double curFit   );
		void initBestCand (bool learnStageFlag );
		void eliminateError ( PopulationT<double> population, unsigned index);
		void sortGroupInfo ();
		unsigned countFailGroupNum ();
		void getPriorInterStage ();
		void combineGroup ( unsigned group1, unsigned group2 );
		void sampleLearnStage ( );
		void binSearchLearnStage (  ); 
		int findInteractPosition ( IndividualT<double> indiv1, IndividualT<double> indiv2, unsigned indexI ); 
		bool testInteraction ( IndividualT<double> indiv1, IndividualT<double> indiv2, unsigned indexI ); 
		void RandomSampleGenDef (); 
		bool TestInterWalk (  IndividualT<double> indiv2, unsigned indexI, IndividualT<double> &localBest , vector< IndividualT<double> > &improveIndivVec); 
		bool TestInterWalk (  IndividualT<double> indiv2, unsigned indexI, IndividualT<double> &localBest); 
		void RandomWalkGenDef (  ); 
		void BinSearchRandWalk (  ); 
		int findInterPosWalk ( IndividualT<double> &localBest, IndividualT<double> indiv2, unsigned indexI , vector< IndividualT<double> > &improveIndivVec); 

	public:
		CCVIL(RunParameter* runParam);
		~CCVIL();
		void run();
		void setObjFunc(Benchmarks* infp){fp = infp;};
};
#endif
