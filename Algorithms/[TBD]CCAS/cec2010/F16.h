#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F16 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F16();
		double compute(double* x);
		double compute(vector<double> x);
		~F16();
	};
}


