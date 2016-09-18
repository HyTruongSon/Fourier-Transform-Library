// Software: Discrete Fourier Transform 2D (for real signals only)
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

// M: number of rows
// N: number of columns
int M, N;

// Normalization constant
double c;

// The original signal
double **s;

// Its Fourier transform
double **Re_F;
double **Im_F;

// Temparary matrix
double **Re_P;
double **Im_P;

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
        std::cerr << "The number of input parameters must be exactly 1 (only for the real signal)!" << std::endl;
        return;
    }
    
    if (nOutputs != 2) {
        std::cerr << "The number of output parameters must be exactly 2 (real and imaginary parts)!" << std::endl;
        return;
    }
    
    M = mxGetM(input_pointers[0]);
    N = mxGetN(input_pointers[0]);
    
    // Memory allocation
    s = new double* [M];
    Re_F = new double* [M];
    Im_F = new double* [M];
    Re_P = new double* [M];
    Im_P = new double* [M];
    
    for (int row = 0; row < M; ++row) {
        s[row] = new double [N];
        Re_F[row] = new double [N];
        Im_F[row] = new double [N];
        Re_P[row] = new double [N];
        Im_P[row] = new double [N];
    }
    
    // Initialization
    vector2matrix(mxGetPr(input_pointers[0]), M, N, s);
    
    // Normalization constant
    c = 1.0 / sqrt(M * N);
    
    // Discrete Fourier Transform - Precomputation
    for (int u = 0; u < M; ++u) {
        for (int y = 0; y < N; ++y) {
            Re_P[u][y] = 0.0;
            Im_P[u][y] = 0.0;
            
            double t = 2.0 * M_PI * (double)(y) / (double)(N);
            for (int v = 0; v < N; ++v) {
                double alpha = t * (double)(v);
                Re_P[u][y] += s[u][v] * cos(alpha);
                Im_P[u][y] -= s[u][v] * sin(alpha);
            }
        }
    }
    
    // Discrete Fourier Transform - Computation
    for (int x = 0; x < M; ++x) {
        for (int y = 0; y < N; ++y) {
            Re_F[x][y] = 0.0;
            Im_F[x][y] = 0.0;
            
            double t = 2.0 * M_PI * (double)(x) / (double)(M);
            for (int u = 0; u < M; ++u) {
                double re = Re_P[u][y];
                double im = Im_P[u][y];
                double alpha = t * (double)(u);
                double cosine = cos(alpha);
                double sine = sin(alpha);
                
                Re_F[x][y] += re * cosine - im * sine;
                Im_F[x][y] += im * cosine + re * sine;
            }
            
            Re_F[x][y] *= c;
            Im_F[x][y] *= c;
        }
    }
    
    // Return outputs
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(Re_F, M, N, mxGetPr(output_pointers[0]));
    
    output_pointers[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(Im_F, M, N, mxGetPr(output_pointers[1]));
}