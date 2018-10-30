#include "F1.h"

/**
 * Shifted Elliptic Function
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

	F1::F1() :Benchmarks(){
		m_havenextGaussian = 0;
		Ovector = NULL;
		minX = -100;
		maxX = 100;
		ID = 1;
		lookup = lookupprepare(dimension);
	}

	F1::~F1(){
		delete[] Ovector;
		delete[] lookup;
	}

	double F1::compute(double* x) {
		double result;
		int    i;

		if (Ovector == NULL) {
			Ovector = createShiftVector(dimension, minX, maxX);
		}

		for (i = dimension - 1; i >= 0; i--) {
			anotherz[i] = x[i] - Ovector[i];
		}

		result = elliptic(anotherz, dimension);
		return(result);
	}


	double F1::compute(vector<double> x){
		double result;
		int    i;

		if (Ovector == NULL) {
			Ovector = createShiftVector(dimension, minX, maxX);
		}

		for (i = dimension - 1; i >= 0; i--) {
			anotherz[i] = x[i] - Ovector[i];
		}

		result = elliptic(anotherz, dimension);
		return(result);
	}
}