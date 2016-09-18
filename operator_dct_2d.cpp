// Software: Discrete Cosine Transform Operator 2D for JPEG
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

double alpha(int u) {
    if (u == 0) {
        return 1.0 / sqrt(2.0);
    }
    return 1.0;
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2 (Number of rows, and Number of columns)!" << std::endl;
        return;
    }
    
    // Number of rows - First input parameter
    int M = mxGetScalar(input_pointers[0]);
    
    // Number of columns - Second input parameter
    int N = mxGetScalar(input_pointers[1]);
    
    // Memory allocation
    /*
    double **block = new double* [M];
    for (int i = 0; i < M; ++i) {
        block[i] = new double [N];
    }
    */
    
    double **result = new double* [M * M];
    for (int i = 0; i < M * M; ++i) {
        result[i] = new double [N * N];
    }
    
    // Computation
    for (int u = 0; u < M; ++u) {
        for (int v = 0; v < N; ++v) {
            double c = 0.25 * alpha(u) * alpha(v);
            double t1 = (double)(u) * M_PI / (2.0 * M);
            double t2 = (double)(v) * M_PI / (2.0 * N);
            
            for (int x = 0; x < M; ++x) {
                for (int y = 0; y < N; ++y) {
                    result[u * M + x][v * N + y] = c * cos((2.0 * x + 1) * t1) * cos((2.0 * y + 1) * t2);
                }
            }
        }
    }
    
    // Return the result
    output_pointers[0] = mxCreateDoubleMatrix(M * M, N * N, mxREAL);
    matrix2vector(result, M * M, N * N, mxGetPr(output_pointers[0]));
}
