#include "F11.h"
#include <stdio.h>

/**
 * D/2m-group Shifted and m-rotated Ackley’s Function
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
	F11::F11() :Benchmarks(){
		m_havenextGaussian = 0;
		Ovector = NULL;
		minX = -32;
		maxX = 32;
		ID = 11;
	}

	F11::~F11(){
		delete[] Ovector;
		delete[] Pvector;
		delete[] RotMatrix;
	}

	double F11::compute(double*x){
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
			result += rot_ackley(anotherz, nonSeparableGroupSize, k);
		}
		//printf("Rot Ackley = %1.16E\n", result);

		double sepSum = ackley(anotherz, dimension, 2);
		//printf("Separable Ackley = %1.16E\n", sepSum);

		result += sepSum;
		return(result);
	}


	double F11::compute(vector<double> x){
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
			result += rot_ackley(anotherz, nonSeparableGroupSize, k);

		}
		//	printf("Rot Ackley = %1.16E\n", result);

		double sepSum = ackley(anotherz, dimension, 2);
		//	printf("Separable Ackley = %1.16E\n", sepSum);

		result += sepSum;
		return(result);
	}
}