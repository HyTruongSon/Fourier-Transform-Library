// Software: Discrete Fourier Transform 1D (for real and complex signals)
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <ctime>

#include "mex.h"

using namespace std;

int N;
double *Re;
double *Im;
double *Amplitude;

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "Exactly 2 input parameters!" << std::endl;
        return;
    }
    
    if (nOutputs != 1) {
        std::cerr << "Exactly 1 output parameters!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) > (size_t)(1)) && (mxGetN(input_pointers[0]) > (size_t)(1))) {
        std::cerr << "Parameter as vector only!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) || (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1]))) {
        std::cerr << "The size of the real part and the imaginary part must be the same!" << std::endl;
        return;
    }
    
    if (mxGetM(input_pointers[0]) > (size_t)(1)) {
        N = mxGetM(input_pointers[0]);
    } else {
        N = mxGetN(input_pointers[1]);
    }
    
    Re = mxGetPr(input_pointers[0]);
    Im = mxGetPr(input_pointers[1]);
    
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Amplitude = mxGetPr(output_pointers[0]);
    
    for (int f = 0; f < N; ++f) {
        Amplitude[f] = sqrt(Re[f] * Re[f] + Im[f] * Im[f]);
    }
}