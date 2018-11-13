#include "F9.h"
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

F9::F9(RunParameter* runParam):Benchmarks(runParam){
	dimension = runParam->dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 9;
	lookup = lookupprepare(nonSeparableGroupSize);
	lookup2 = lookupprepare(dimension/2);
}

F9::F9():Benchmarks(){
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 9;
	lookup = lookupprepare(nonSeparableGroupSize);
	lookup2 = lookupprepare(dimension/2);
}

F9::~F9(){
	delete[] Ovector;
	delete[] Pvector;
	delete[] lookup;
	delete[] lookup2;
	// delete 2D array
	int i;
	for(i=0;i<dimension/(2*nonSeparableGroupSize);i++){
		delete[] MultiRotMatrix1D[i];
	}
	delete[] MultiRotMatrix1D;
}

double F9::compute(double*x){
	int k, i;
	double result=0.0;

	if(Ovector==NULL){
		Ovector=createShiftVector(dimension,minX,maxX);
		Pvector=createPermVector(dimension);
		MultiRotMatrix1D=createMultiRotateMatrix1D(nonSeparableGroupSize,dimension/(2*nonSeparableGroupSize));

		/* 
		 * print the multi rotated matrix 
		printf("\n\n\n print the multi rotated matrix\n\n\n");
		for (k = 0; k<dimension/(2*nonSeparableGroupSize); k++){
		printf("\n matrix %d: \n", k+1);
			for (i = 0; i<nonSeparableGroupSize*nonSeparableGroupSize; i++){
				printf("%1.20E\t", MultiRotMatrix1D[k][i]);
			}
		}
		 */
	}

	for( i=0;i<dimension;i++){
		anotherz[i]=x[i]-Ovector[i];
	}

	//
	//	printf ( "Pvector\n" );
	//	for(i=0;i<dimension;i++){
	//		printf ( "%d\n", Pvector[i] );
	//	}

	for(k=1;k<=dimension/(2*nonSeparableGroupSize);k++){
		result+=rot_elliptic(anotherz,nonSeparableGroupSize,k);
	}

//	printf("Rotated Part = %1.20E\n", result);
//	printf("Non-Rotated Part = %1.20E\n", elliptic(anotherz, dimension, 2));

	result+=elliptic(anotherz, dimension, 2);

	return(result);
}

double F9::compute(vector<double> x){
	int i,k;
	double result=0.0;

	if(Ovector==NULL){
		Ovector=createShiftVector(dimension,minX,maxX);
		Pvector=createPermVector(dimension);
		MultiRotMatrix1D=createMultiRotateMatrix1D(nonSeparableGroupSize,dimension/(2*nonSeparableGroupSize));

		/* 
		 * print the multi rotated matrix 
		printf("\n\n\n print the multi rotated matrix\n\n\n");
		for (k = 0; k<dimension/(2*nonSeparableGroupSize); k++){
		printf("\n matrix %d: \n", k+1);
			for (i = 0; i<nonSeparableGroupSize*nonSeparableGroupSize; i++){
				printf("%1.20E\t", MultiRotMatrix1D[k][i]);
			}
		}
		 */
	}

	for(i=0;i<dimension;i++){
		anotherz[i]=x[i]-Ovector[i];
	}

	for(k=1;k<=dimension/(2*nonSeparableGroupSize);k++){
		result+=rot_elliptic(anotherz,nonSeparableGroupSize,k);
	}

//	printf("Rotated Part = %1.20E\n", result);

	result+=elliptic(anotherz, dimension, 2);
	return(result);
}
