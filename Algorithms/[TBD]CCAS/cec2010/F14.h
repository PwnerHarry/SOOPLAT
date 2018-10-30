#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F14 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F14();
		double compute(double* x);
		double compute(vector<double> x);
		~F14();
	};
}


