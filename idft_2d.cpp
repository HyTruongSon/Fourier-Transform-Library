// Software: Inverse Discrete Fourier Transform 2D (for real original signals only)
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
double **Re_H;
double **Im_H;

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
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2 (real and imaginary parts)!" << std::endl;
        return;
    }
    
    if (nOutputs != 1) {
        std::cerr << "The number of output parameters must be exactly 1 (only for the real signal)!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) || (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1]))) {
        std::cerr << "The real and imaginary parts must have the same size!" << std::endl;
        return;
    }
    
    M = mxGetM(input_pointers[0]);
    N = mxGetN(input_pointers[0]);
    
    // Memory allocation
    s = new double* [M];
    Re_F = new double* [M];
    Im_F = new double* [M];
    Re_H = new double* [M];
    Im_H = new double* [M];
    
    for (int row = 0; row < M; ++row) {
        s[row] = new double [N];
        Re_F[row] = new double [N];
        Im_F[row] = new double [N];
        Re_H[row] = new double [N];
        Im_H[row] = new double [N];
    }
    
    // Initialization
    vector2matrix(mxGetPr(input_pointers[0]), M, N, Re_F);
    vector2matrix(mxGetPr(input_pointers[1]), M, N, Im_F);
    
    // Normalization constant
    c = 1.0 / sqrt(M * N);
    
    // Discrete Fourier Transform - Precomputation
    for (int x = 0; x < M; ++x) {
        for (int v = 0; v < N; ++v) {
            Re_H[x][v] = 0.0;
            Im_H[x][v] = 0.0;
            
            double t = 2.0 * M_PI * (double)(v) / (double)(N);
            for (int y = 0; y < N; ++y) {
                double alpha = t * (double)(y);
                double cosine = cos(alpha);
                double sine = sin(alpha);
                double re = Re_F[x][y];
                double im = Im_F[x][y];
            
                Re_H[x][v] += re * cosine - im * sine;
                Im_H[x][v] += im * cosine + re * sine;
            }
        }
    }
    
    // Discrete Fourier Transform - Computation
    for (int u = 0; u < M; ++u) {
        for (int v = 0; v < N; ++v) {
            s[u][v] = 0.0;
            
            double t = 2.0 * M_PI * (double)(u) / (double)(M);
            for (int x = 0; x < M; ++x) {
                double alpha = t * (double)(x);
                s[u][v] += Re_H[x][v] * cos(alpha) - Im_H[x][v] * sin(alpha);
            }
            
            s[u][v] *= c;
        }
    }
    
    // Return outputs
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(s, M, N, mxGetPr(output_pointers[0]));
}