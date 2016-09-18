// Software: Normal polynomial multiplication
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time Complexity: O(N^2)
// Space Complexity: O(N)

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

int nA, nB, nC;
double *A, *B, *C;

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2 (polynomials)!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) > size_t(1)) && (mxGetN(input_pointers[0]) > size_t(1))) {
        std::cerr << "The input parameters must be vectors!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[1]) > size_t(1)) && (mxGetN(input_pointers[1]) > size_t(1))) {
        std::cerr << "The input parameters must be vectors!" << std::endl;
        return;
    }
    
    A = mxGetPr(input_pointers[0]);
    if (mxGetM(input_pointers[0]) > size_t(1)) {
        nA = mxGetM(input_pointers[0]);
    } else {
        nA = mxGetN(input_pointers[0]);
    }
    
    B = mxGetPr(input_pointers[1]);
    if (mxGetM(input_pointers[1]) > size_t(1)) {
        nB = mxGetM(input_pointers[1]);
    } else {
        nB = mxGetN(input_pointers[1]);
    }
    
    nC = nA + nB - 1;
    output_pointers[0] = mxCreateDoubleMatrix(nC, 1, mxREAL);
    C = mxGetPr(output_pointers[0]);
    
    for (int i = 0; i < nC; ++i) {
        C[i] = 0.0;
        int L = max(0, i - nB + 1);
        int R = min(i, nA - 1);
        for (int j = L; j <= R; ++j) {
            C[i] += A[j] * B[i - j];
        }
    }
}