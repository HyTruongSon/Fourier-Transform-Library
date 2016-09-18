// Software: Inverse Discrete Fourier Transform Operator 1D (use for further computation)
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

void vector2matrix(double *input, int nRows, int nCols, double **output) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            output[i][j] = input[j * nRows + i];
        }
    }
}

void matrix2vector(double **input, int nRows, int nCols, double *output) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            output[j * nRows + i] = input[i][j];
        }
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 1) {
        std::cerr << "The number of input parameters must be exactly 1!" << std::endl;
        return;
    }
    
    // Input the dimension
    int N = mxGetScalar(input_pointers[0]);
        
    // Memory allocation
    double **Re_w = new double* [N];
    double **Im_w = new double* [N];
    
    for (int u = 0; u < N; ++u) {
        Re_w[u] = new double [N];
        Im_w[u] = new double [N];
    }
    
    // Computation
    for (int u = 0; u < N; ++u) {
        double t = 2.0 * M_PI * (double)(u) / (double)(N);
        for (int v = u; v < N; ++v) {
            double alpha = t * (double)(v);
            
            Re_w[u][v] = cos(alpha);
            Im_w[u][v] = sin(alpha);
            
            Re_w[v][u] = Re_w[u][v];
            Im_w[v][u] = Im_w[u][v];
        }
    }
    
    // Return results
    output_pointers[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    matrix2vector(Re_w, N, N, mxGetPr(output_pointers[0]));
    
    output_pointers[1] = mxCreateDoubleMatrix(N, N, mxREAL);
    matrix2vector(Im_w, N, N, mxGetPr(output_pointers[1]));
}