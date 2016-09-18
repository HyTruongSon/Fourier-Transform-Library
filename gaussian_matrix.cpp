// Software: Gaussian matrix
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
        std::cerr << "The number of input parameters must be 1 (a vector of the radius and variance)!" << std::endl;
        return;
    }
    
    double *parameters = mxGetPr(input_pointers[0]);
    int r = (int)(parameters[0]);
    double o2 = parameters[1];
    
    if ((r <= 0) || (o2 <= 0.0)) {
        std::cerr << "The radius and the variance must be greater than 0!" << std::endl;
        return;
    }
    
    int N = 2 * r + 1;
    double **matrix = new double* [N];
    for (int i = 0; i < N; ++i) {
        matrix[i] = new double [N];
    }
    
    double c = 1.0 / (2.0 * M_PI * o2);
    for (int x = 0; x <= r; ++x) {
        for (int y = 0; y <= r; ++y) {
            matrix[r + x][r + y] = c * exp(- (x * x + y * y) / (2.0 * o2));
            matrix[r - x][r + y] = matrix[r + x][r + y];
            matrix[r - x][r - y] = matrix[r + x][r + y];
            matrix[r + x][r - y] = matrix[r + x][r + y];
        }
    }
    
    output_pointers[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    matrix2vector(matrix, N, N, mxGetPr(output_pointers[0]));
}