// Software: Convert from vector form to matrix form
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
    if (nInputs != 3) {
        std::cerr << "The number of input parameters must be exactly 3 (the vector, the number of rows, and the number of columns)!" << std::endl;
        return;
    }
    
    // Dimensions
    int nRows = mxGetScalar(input_pointers[1]);
    int nCols = mxGetScalar(input_pointers[2]);

    // Return the output
    output_pointers[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    copy(mxGetPr(input_pointers[0]), nRows * nCols, mxGetPr(output_pointers[0]));
}