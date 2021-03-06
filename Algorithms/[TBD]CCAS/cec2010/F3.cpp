#include "F3.h"

/**
 * Shifted Ackley's Function
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
	F3::F3() :Benchmarks(){
		m_havenextGaussian = 0;
		Ovector = NULL;
		minX = -32;
		maxX = 32;
		ID = 3;
	}

	F3::~F3(){
		delete[] Ovector;
	}


	double F3::compute(double*x){
		int    i;
		double result;

		if (Ovector == NULL) {
			Ovector = createShiftVector(dimension, minX, maxX);
		}

		for (i = dimension - 1; i >= 0; i--) {
			anotherz[i] = x[i] - Ovector[i];
		}

		result = ackley(anotherz, dimension);
		return(result);
	}


	double F3::compute(vector<double> x){
		int    i;
		double result;

		if (Ovector == NULL) {
			Ovector = createShiftVector(dimension, minX, maxX);
		}

		for (i = dimension - 1; i >= 0; i--) {
			anotherz[i] = x[i] - Ovector[i];
		}

		result = ackley(anotherz, dimension);
		return(result);
	}
}