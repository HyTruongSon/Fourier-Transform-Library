// Software: Convert from matrix form to vector form
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

void copy(double *input, int N, double *output) {
    for (int i = 0; i < N; ++i) {
        output[i] = input[i];
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 1) {
        std::cerr << "The number of input parameters must be exactly 1!" << std::endl;
        return;
    }
    
    // Dimensions
    int N = mxGetM(input_pointers[0]) * mxGetN(input_pointers[0]);

    // Return the output
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    copy(mxGetPr(input_pointers[0]), N, mxGetPr(output_pointers[0]));
}