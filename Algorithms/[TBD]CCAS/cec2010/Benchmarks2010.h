#pragma once

#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdint>
#include "Fitness.h"

using namespace std;


namespace CEC2010 {

	struct IndexMap{
		unsigned arrIndex1;
		unsigned arrIndex2;
	};

	class Benchmarks : public Fitness
	{
	protected:
		int next(int bits);
		int nextInt(int n);
		double nextDouble();
		double nextGaussian();
	
		double* createShiftVector(int dim, double min, double max);
		int* createPermVector(int dim);
		double** createRotMatrix(int dim);
		double* createRotMatrix1D(int dim);
		double** createMultiRotateMatrix1D(int dim, int num);

		double* lookupprepare(int dim);

		// Basic mathematical functions' declaration
		double* multiply(double*vector, double*matrix, int dim);
		double elliptic(double*x, int dim);
		double elliptic(double*x, int dim, int k);
		double rastrigin(double*x, int dim);
		double rastrigin(double *x, int dim, int k);
		double ackley(double*x, int dim);
		double ackley(double*x, int dim, int k);
		double rot_elliptic(double*x, int dim);
		double rot_elliptic(double*x, int dim, int k);
		double rot_rastrigin(double*x, int dim);
		double rot_rastrigin(double *x, int dim, int k);
		double rot_ackley(double*x, int dim);
		double rot_ackley(double*x, int dim, int k);
		double schwefel(double*x, int dim);
		double schwefel(double*x, int dim, int k);
		double sphere(double*x, int dim);
		double sphere(double*x, int dim, int k);
		double rosenbrock(double*x, int dim);
		double rosenbrock(double*x, int dim, int k);
		unsigned convertMatrixToArrayIndex(unsigned i, unsigned j);
		void createIndexMapping();

		int64_t M;
		int64_t A;
		int64_t m_seed;
		int64_t MASK;
		double m_nextGaussian;
		bool  m_havenextGaussian;
		bool setOvectorToZero;

		double *Ovector;
		int*    Pvector;
		double* RotMatrix;
		double** MultiRotMatrix1D;
		double *lookup;
		double *lookup2;

		double* anotherz;
		double* anotherz1;
		double* anotherz2;

		vector<bool> interArray;

		// running time setting for benchmarks
		int minX;
		int maxX;
		int dimension;
		int nonSeparableGroupSize;
		int64_t functionInitRandomSeed;
		struct IndexMap *indexMap;
		unsigned arrSize;

	public:
		Benchmarks();
		virtual ~Benchmarks();
		virtual void createIdealGroups();
		void setMinX(int);
		void setMaxX(int);
		void setSeed(int64_t);
		void setDimension(int);
		void setNonSeparableGroupSize(int);
		vector<bool> getInterArray();
		void ArrToMat(unsigned I1, unsigned I2, unsigned &matIndex);
		void MatToArr(unsigned &I1, unsigned &I2, unsigned matIndex);
		virtual double getMinX();
    virtual double getMaxX();
		virtual unsigned getID();
    virtual unsigned  getDimension();
	};

}
