#include "F14.h"
#include <stdio.h>

/**
 * Single-group Shifted and m-rotated Elliptic Function
 *
 * as defined in "Benchmark Functions for the CEC'2010 Special Session
 * and Competition on Large-Scale Global Optimization" by Ke Tang,
 * Xiaodong Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise
 * published as technical report on January 8, 2010 at Nature Inspired
 * Computation and Applications Laboratory (NICAL), School of Computer
 * Science and Technology, University of Science and Technology of China,
 * Hefei, Anhui, China.
 */

F14::F14(RunParameter* runParam):Benchmarks(runParam){
	dimension = runParam->dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 14;

	lookup = lookupprepare(nonSeparableGroupSize);

	Ovector=createShiftVector(dimension,minX,maxX);

//	Pvector = (int*)malloc(sizeof(int) * dimension);
//	for (unsigned i = 0; i<(unsigned)dimension; i++){
//		Pvector[i] = i;	
//	}


	Pvector=createPermVector(dimension);

	MultiRotMatrix1D = createMultiRotateMatrix1D(nonSeparableGroupSize,dimension/(nonSeparableGroupSize));

	generateInterArray ( );
}

F14::F14():Benchmarks(){
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 14;
	lookup = lookupprepare(nonSeparableGroupSize);

	Ovector=createShiftVector(dimension,minX,maxX);
	Pvector=createPermVector(dimension);
	MultiRotMatrix1D=createMultiRotateMatrix1D(nonSeparableGroupSize,dimension/(nonSeparableGroupSize));
}

F14::~F14(){
	delete[] Ovector;
	delete[] Pvector;
	delete[] lookup;
	// delete 2D array
	int i;
	for(i=0;i<dimension/(nonSeparableGroupSize);i++){
		delete[] MultiRotMatrix1D[i];
	}
	delete[] MultiRotMatrix1D;
}

double F14::compute(double*x){
	int i,k;
	double result=0.0;

	for(i=0;i<dimension;i++){
		anotherz[i]=x[i]-Ovector[i];
	}

	for(k=1;k<=dimension/(nonSeparableGroupSize);k++){
		result+=rot_elliptic(anotherz,nonSeparableGroupSize,k);
	}

	return(result);
}

double F14::compute(vector<double> x){
	int i,k;
	double result=0.0;

	for(i=0;i<dimension;i++){
		anotherz[i]=x[i]-Ovector[i];
	}

	for(k=1;k<=dimension/(nonSeparableGroupSize);k++){
		result+=rot_elliptic(anotherz,nonSeparableGroupSize,k);
	}

	return(result);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  F14::generateInterArray
 *  Description:  
 * =====================================================================================
 */
	void
F14::generateInterArray ( )
{
	// initialize the basic structure
	for (unsigned i=0; i<(unsigned)dimension*(dimension-1)/2; i++){
		interArray.push_back(false);
	}

//	printf ( "Print P vector\n" );
//	for (unsigned i=0; i<(unsigned)dimension; i++){
//		printf ( "%d\t", Pvector[i] );
//	}

	// assign values
	unsigned baseIndex=0, compIndex=0;
	for (unsigned i=0; i<(unsigned)dimension/nonSeparableGroupSize; i++){
		for (unsigned j=0; j<(unsigned)nonSeparableGroupSize; j++){
			baseIndex =	Pvector[i*nonSeparableGroupSize+j];
			for (unsigned k=j+1; k<(unsigned)nonSeparableGroupSize; k++){
				compIndex = Pvector[i*nonSeparableGroupSize+k];
				if (baseIndex < compIndex){
//					printf ( "Mat: smallIndex %d, bigIndex %d; Arr: %d\n", baseIndex, compIndex, convertMatrixToArrayIndex(baseIndex, compIndex));
					interArray[convertMatrixToArrayIndex(baseIndex, compIndex)] = true;
				}else{
//					printf ( "Mat: smallIndex %d, bigIndex %d; Arr: %d\n", compIndex, baseIndex, convertMatrixToArrayIndex(compIndex, baseIndex));
//					printf ( "%d\n", convertMatrixToArrayIndex(compIndex, baseIndex));
					interArray[convertMatrixToArrayIndex( compIndex, baseIndex)] = true;
				}
			}
		}
	}
}		/* -----  end of function F14::generateMat  ----- */
