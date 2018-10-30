#ifndef _HEADER_H
#define _HEADER_H


#include <iostream>
#include <cstring>
#include <cstdio>
#include <EALib/PopulationT.h>

// optimization algorithms
#include "CCVIL.h"

// benchmark set header files
#include "benchmark/F1.h"
#include "benchmark/F2.h"
#include "benchmark/F3.h"
#include "benchmark/F4.h"
#include "benchmark/F5.h"
#include "benchmark/F6.h"
#include "benchmark/F7.h"
#include "benchmark/F8.h"
#include "benchmark/F9.h"
#include "benchmark/F10.h"
#include "benchmark/F11.h"
#include "benchmark/F12.h"
#include "benchmark/F13.h"
#include "benchmark/F14.h"
#include "benchmark/F15.h"
#include "benchmark/F16.h"
#include "benchmark/F17.h"
#include "benchmark/F18.h"
#include "benchmark/F19.h"
#include "benchmark/F20.h"


using namespace std;

Benchmarks* generateFuncObj(RunParameter* runParam, int funcID);
Benchmarks* generateFuncObj(int funcID);

#endif
