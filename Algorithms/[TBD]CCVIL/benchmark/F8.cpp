#include "F8.h"
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

F8::F8(RunParameter* runParam):Benchmarks(runParam){
	dimension = runParam->dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 8;
}

F8::F8():Benchmarks(){
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 8;
}

F8::~F8(){
	delete[] Ovector;
	delete[] Pvector;
}

double F8::compute(double* x){
	int    m = nonSeparableGroupSize;
	int    i;
	double result;

	if(Ovector == NULL) {
		Ovector = createShiftVector(dimension,minX,maxX-1);

//		printf("\n\n\nO vector\n\n\n");
//		for (i = 0; i<dimension; i++){
//			printf("%f\t",Ovector[i]);
//		}

		Pvector = createPermVector(dimension);
		
//TODO: Neeed to change back to random one ****************************************************************************
//		Pvector = (int*)malloc(sizeof(int) * dimension);
//		for (i = 0; i<dimension; i++){
//			Pvector[i] = i;	
//		}

		/*
		printf("\n\n\nP vector\n\n\n");
		for (i = 0; i<dimension; i++){
			printf("%d\t",Pvector[i]);
		}
		*/
	}

	for(i = 0; i < dimension; i++) {
		anotherz[i] = x[i] - Ovector[i];
	}

	for(i = 0; i < m; i++) {
		anotherz1[i] = anotherz[Pvector[i]];
	}

	for(i = m; i < dimension; i++) {
		anotherz2[i - m] = anotherz[Pvector[i]];
	}

//	printf("\n\n\nanotherz1\n\n\n");
//	for (i = 0; i<m; i++){
//		printf("%f\t",anotherz1[i]);
//	}

	result = rosenbrock(anotherz1,m) * 1e6 + sphere(anotherz2,dimension - m);

//	printf("Rosenbrock Part = %1.16E\n", rosenbrock(anotherz1,m) * 1e6);
//	printf("Sphere Part = %1.16E\n", sphere(anotherz2,dimension - m));

	return(result);
}


double F8::compute(vector<double> x){
	int    m = nonSeparableGroupSize;
	int    i;
	double result;

	if(Ovector == NULL) {

		Ovector = createShiftVector(dimension,minX,maxX-1);

		//		printf("\n\n\nO vector\n\n\n");
		//		for (i = 0; i<dimension; i++){
		//			printf("%f\t",Ovector[i]);
		//		}

		Pvector = createPermVector(dimension);
		
		//		//TODO: Neeed to change back to random one ****************************************************************************
		//		Pvector = (int*)malloc(sizeof(int) * dimension);
		//		for (i = 0; i<dimension; i++){
		//			Pvector[i] = i;	
		//		}

		//		printf("\n\n\nP vector\n\n\n");
		//		for (i = 0; i<dimension; i++){
		//			printf("%d\t",Pvector[i]);
		//		}
		//		printf ( "\n" );
		//
		//		printf ( "dimension = %d\n", dimension );
		//		printf ( "m = %d\n", m );
	}

	for(i = 0; i < dimension; i++) {
		anotherz[i] = x[i] - Ovector[i];
	}

	for(i = 0; i < m; i++) {
		anotherz1[i] = anotherz[Pvector[i]];
	}

	for(i = m; i < dimension; i++) {
		anotherz2[i - m] = anotherz[Pvector[i]];
	}
			
	result = rosenbrock(anotherz1,m) * 1e6 + sphere(anotherz2,dimension - m);

	//	printf("Rosenbrock Part = %1.16E\n", rosenbrock(anotherz1,m) * 1e6);
	//	printf("Sphere Part = %1.16E\n", sphere(anotherz2,dimension - m));

	return(result);
}
