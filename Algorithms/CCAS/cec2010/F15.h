#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F15 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F15();
		double compute(double* x);
		double compute(vector<double> x);
		~F15();
	};
}

