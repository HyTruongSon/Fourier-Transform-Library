// Software: Discrete Fourier Transform Opeartor 1D (only use for visualization)
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
    
    // Computation - We need to shift such that (0, 0) is in the middle
    int half = N / 2;
    for (int u = 0; u < N; ++u) {
        int i = (u + half) % N;
        double t = 2.0 * M_PI * (double)(u) / (double)(N);
        for (int v = u; v < N; ++v) {
            int j = (v + half) % N;
            double alpha = t * (double)(v);
            
            Re_w[i][j] = cos(alpha);
            Im_w[i][j] = - sin(alpha);
            
            Re_w[j][i] = Re_w[i][j];
            Im_w[j][i] = Im_w[i][j];
        }
    }
    
    // Return results
    output_pointers[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    matrix2vector(Re_w, N, N, mxGetPr(output_pointers[0]));
    
    output_pointers[1] = mxCreateDoubleMatrix(N, N, mxREAL);
    matrix2vector(Im_w, N, N, mxGetPr(output_pointers[1]));
}