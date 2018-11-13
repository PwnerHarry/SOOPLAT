#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F19 :public Benchmarks{
	protected:
	public:
		F19();
		double compute(double* x);
		double compute(vector<double> x);
		~F19();
	};
}



