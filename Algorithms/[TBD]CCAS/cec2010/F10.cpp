#include "F10.h"
#include <stdio.h>

/**
 * Single-group Shifted and m-rotated Rastrigin’s Function
 *
 * as defined in "Benchmark Functions for the CEC'2010 Special Session
 * and Competition on Large-Scale Global Optimization" by Ke Tang,
 * Xiaodong Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise
 * published as technical report on January 8, 2010 at Nature Inspired
 * Computation and Applications Laboratory (NICAL), School of Computer
 * Science and Technology, University of Science and Technology of China,
 * Hefei, Anhui, China.
 */
namespace CEC2010 {
	F10::F10() :Benchmarks(){
		m_havenextGaussian = 0;
		Ovector = NULL;
		minX = -5;
		maxX = 5;
		ID = 10;
	}

	F10::~F10(){
		delete[] Ovector;
		delete[] Pvector;
		delete[] RotMatrix;
	}

	double F10::compute(double*x){
		int i, k;
		double result = 0.0;

		if (Ovector == NULL)
		{
			Ovector = createShiftVector(dimension, minX, maxX);
			Pvector = createPermVector(dimension);
			RotMatrix = createRotMatrix1D(nonSeparableGroupSize);
		}
		for (i = 0; i < dimension; i++)
		{
			anotherz[i] = x[i] - Ovector[i];
		}
		for (k = 1; k <= dimension / (2 * nonSeparableGroupSize); k++)
		{
			result += rot_rastrigin(anotherz, nonSeparableGroupSize, k);
		}

		//	printf("Rot Rastrigin = %1.16E\n", result);

		result += rastrigin(anotherz, dimension, 2);
		return(result);
	}

	double F10::compute(vector<double> x){
		int i, k;
		double result = 0.0;

		if (Ovector == NULL)
		{
			Ovector = createShiftVector(dimension, minX, maxX);
			Pvector = createPermVector(dimension);
			RotMatrix = createRotMatrix1D(nonSeparableGroupSize);
		}
		for (i = 0; i < dimension; i++)
		{
			anotherz[i] = x[i] - Ovector[i];
		}
		for (k = 1; k <= dimension / (2 * nonSeparableGroupSize); k++)
		{
			result += rot_rastrigin(anotherz, nonSeparableGroupSize, k);
		}

		//	printf("Rot Rastrigin = %1.16E\n", result);

		result += rastrigin(anotherz, dimension, 2);
		return(result);
	}
}