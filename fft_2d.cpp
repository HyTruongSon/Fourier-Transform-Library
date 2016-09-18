// Software: Fast Fourier Transform 2D (for real signals only)
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time complexity: O(N^2logN)
// Space complexity: O(N^2)

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

// Fast Fourier Transform 1D for real signals
void FFT_real(double *Signal, double *Re_F, double *Im_F, int N, int t) {
    if (N == 1) {
        Re_F[0] = Signal[0];
        Im_F[0] = 0.0;
        return;
    }
    
    int half = N / 2;
    FFT_real(Signal, Re_F, Im_F, half, 2 * t);
    FFT_real(Signal + t, Re_F + half, Im_F + half, half, 2 * t);
    
    for (int k = 0; k < half; ++k) {
        double r1 = Re_F[k];
        double i1 = Im_F[k];
        
        double r2 = Re_F[k + half];
        double i2 = Im_F[k + half];
        
        double alpha = - 2.0 * M_PI * (double)(k) / (double)(N);
        double r3 = cos(alpha);
        double i3 = sin(alpha);
        
        double r4 = r2 * r3 - i2 * i3;
        double i4 = r3 * i2 + r2 * i3;
        
        Re_F[k] = r1 + r4;
        Im_F[k] = i1 + i4;
        
        Re_F[k + half] = r1 - r4;
        Im_F[k + half] = i1 - i4;
    }
}

// Fast Fourier Transform 1D for complex signals
void FFT_complex(double *Re_Signal, double *Im_Signal, double *Re_F, double *Im_F, int N, int t) {
    if (N == 1) {
        Re_F[0] = Re_Signal[0];
        Im_F[0] = Im_Signal[0];
        return;
    }
    
    int half = N / 2;
    FFT_complex(Re_Signal, Im_Signal, Re_F, Im_F, half, 2 * t);
    FFT_complex(Re_Signal + t, Im_Signal + t, Re_F + half, Im_F + half, half, 2 * t);
    
    for (int k = 0; k < half; ++k) {
        double r1 = Re_F[k];
        double i1 = Im_F[k];
        
        double r2 = Re_F[k + half];
        double i2 = Im_F[k + half];
        
        double alpha = - 2.0 * M_PI * (double)(k) / (double)(N);
        double r3 = cos(alpha);
        double i3 = sin(alpha);
        
        double r4 = r2 * r3 - i2 * i3;
        double i4 = r3 * i2 + r2 * i3;
        
        Re_F[k] = r1 + r4;
        Im_F[k] = i1 + i4;
        
        Re_F[k + half] = r1 - r4;
        Im_F[k + half] = i1 - i4;
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
    
    // M: number of rows
    // N: number of columns
    int M = mxGetM(input_pointers[0]);
    int N = mxGetN(input_pointers[0]);
    
    // Normalization constant
    double c = 1.0 / sqrt(M * N);

    // Memory allocation
    double **s = new double* [M];    // The original signal
    double **Re_F = new double* [M]; // Its Fourier transform
    double **Im_F = new double* [M]; // Its Fourier transform
    double **Re_P = new double* [M]; // Temparary matrix
    double **Im_P = new double* [M]; // Temparary matrix
    
    for (int row = 0; row < M; ++row) {
        s[row] = new double [N];
        Re_F[row] = new double [N];
        Im_F[row] = new double [N];
        Re_P[row] = new double [N];
        Im_P[row] = new double [N];
    }
    
    // Initialization
    vector2matrix(mxGetPr(input_pointers[0]), M, N, s);
    
    // Fast Fourier Transform 2D - Precomputation
    /*
    double *A = new double [N];
    double *B = new double [N];
    double *C = new double [N];
    */
     
    for (int u = 0; u < M; ++u) {
        FFT_real(s[u], Re_P[u], Im_P[u], N, 1);
        
        // Alternative
        /*
        for (int v = 0; v < N; ++v) {
            A[v] = s[u][v];
        }
        
        FFT_real(A, B, C, N);
        
        for (int y = 0; y < N; ++y) {
            Re_P[u][y] = B[y];
            Im_P[u][y] = C[y];
        }
        */
    }
    
    // Fast Fourier Transform 2D - Computation
    double *Re_p = new double [M];
    double *Im_p = new double [M];
    double *Re_f = new double [M];
    double *Im_f = new double [M];
     
    for (int y = 0; y < N; ++y) {
        for (int u = 0; u < M; ++u) {
            Re_p[u] = Re_P[u][y];
            Im_p[u] = Im_P[u][y];
        }
        
        FFT_complex(Re_p, Im_p, Re_f, Im_f, M, 1);
        
        for (int x = 0; x < M; ++x) {
            Re_F[x][y] = Re_f[x] * c;
            Im_F[x][y] = Im_f[x] * c;
        }
    }
    
    // Return outputs
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(Re_F, M, N, mxGetPr(output_pointers[0]));
    
    output_pointers[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(Im_F, M, N, mxGetPr(output_pointers[1]));
}