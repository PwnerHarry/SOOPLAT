#include "F20.h"
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

F20::F20(RunParameter* runParam):Benchmarks(runParam){
	dimension = runParam->dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 20;
}

F20::F20():Benchmarks(){
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -100;
	maxX = 100;
	ID = 20;
}

F20::~F20(){
	delete[] Ovector;
}

double F20::compute(double*x){
	int i;
	double result=0.0;

	if(Ovector==NULL)
	{
		Ovector=createShiftVector(dimension,minX,maxX-1);
	}

	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}

	result=rosenbrock(anotherz, dimension);
	return(result);
}

double F20::compute(vector<double> x){
	int i;
	double result=0.0;

	if(Ovector==NULL)
	{
		Ovector=createShiftVector(dimension,minX,maxX-1);
	}

	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}

	result=rosenbrock(anotherz, dimension);
	return(result);
}
