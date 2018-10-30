#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F17 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F17();
		double compute(double* x);
		double compute(vector<double> x);
		~F17();
	};
}


