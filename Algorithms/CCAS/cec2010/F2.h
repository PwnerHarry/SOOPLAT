#pragma once

#include "Benchmarks2010.h"

namespace CEC2010 {
	class F2 :public Benchmarks{
	protected:

	public:
		F2();
		double compute(double* x);
		double compute(vector<double> x);
		~F2();
	};
}

