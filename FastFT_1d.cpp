// Software: Fast Fourier Transform 1D (for real and complex signals)
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time complexity: O(NlogN)
// Space complexity: O(NlogN). Better implementation: fft_1d.cpp - O(N).

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

// type = 1: The original signal is real
// type = 2: The original signal is complex
int type;

// Number of time samples
int N;

// Normalization constant c = 1/sqrt(N)
double c; 

// Fast Fourier Transform 1D for real signals
void FFT_real(double *Signal, double *Re_F, double *Im_F, int N) {
    if (N == 1) {
        Re_F[0] = Signal[0];
        Im_F[0] = 0.0;
        return;
    }
    
    int M = N / 2;
    double *Even = new double [M];
    double *Odd = new double [M];
    
    for (int i = 0; i < M; ++i) {
        Even[i] = Signal[2 * i];
        Odd[i] = Signal[2 * i + 1];
    }
    
    double *Re_F_Even = new double [M];
    double *Im_F_Even = new double [M];
    
    double *Re_F_Odd = new double [M];
    double *Im_F_Odd = new double [M];
    
    FFT_real(Even, Re_F_Even, Im_F_Even, M);
    FFT_real(Odd, Re_F_Odd, Im_F_Odd, M);
    
    for (int k = 0; k < M; ++k) {
        double r1 = Re_F_Even[k];
        double i1 = Im_F_Even[k];
        
        double r2 = Re_F_Odd[k];
        double i2 = Im_F_Odd[k];
        
        double alpha = - 2.0 * M_PI * (double)(k) / (double)(N);
        double r3 = cos(alpha);
        double i3 = sin(alpha);
        
        double r4 = r2 * r3 - i2 * i3;
        double i4 = r3 * i2 + r2 * i3;
        
        Re_F[k] = r1 + r4;
        Im_F[k] = i1 + i4;
        
        Re_F[k + M] = r1 - r4;
        Im_F[k + M] = i1 - i4;
    }
}

// Fast Fourier Transform 1D for complex signals
void FFT_complex(double *Re_Signal, double *Im_Signal, double *Re_F, double *Im_F, int N) {
	if (N == 1) {
        Re_F[0] = Re_Signal[0];
        Im_F[0] = Im_Signal[0];
        return;
    }
    
    int M = N / 2;
    
    double *Re_Even = new double [M];
    double *Im_Even = new double [M];
    
    double *Re_Odd = new double [M];
    double *Im_Odd = new double [M];
    
    for (int i = 0; i < M; ++i) {
        Re_Even[i] = Re_Signal[2 * i];
        Im_Even[i] = Im_Signal[2 * i];
        
        Re_Odd[i] = Re_Signal[2 * i + 1];
        Im_Odd[i] = Im_Signal[2 * i + 1];
    }
    
    double *Re_F_Even = new double [M];
    double *Im_F_Even = new double [M];
    
    double *Re_F_Odd = new double [M];
    double *Im_F_Odd = new double [M];
    
    FFT_complex(Re_Even, Im_Even, Re_F_Even, Im_F_Even, M);
    FFT_complex(Re_Odd, Im_Odd, Re_F_Odd, Im_F_Odd, M);
    
    for (int k = 0; k < M; ++k) {
        double r1 = Re_F_Even[k];
        double i1 = Im_F_Even[k];
        
        double r2 = Re_F_Odd[k];
        double i2 = Im_F_Odd[k];
        
        double alpha = - 2.0 * M_PI * (double)(k) / (double)(N);
        double r3 = cos(alpha);
        double i3 = sin(alpha);
        
        double r4 = r2 * r3 - i2 * i3;
        double i4 = r3 * i2 + r2 * i3;
        
        Re_F[k] = r1 + r4;
        Im_F[k] = i1 + i4;
        
        Re_F[k + M] = r1 - r4;
        Im_F[k + M] = i1 - i4;
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs == 0) {
        std::cerr << "Not enough input parameters!" << std::endl;
        return;
    }
    
    if (nInputs > 2) {
        std::cerr << "Maximum 2 input parameters only!" << std::endl;
        return;
    }
    
    if (nOutputs != 2) {
        std::cerr << "The number of output parameters must be 2 (the real and imaginary parts)!" << std::endl;
        return;
    }
    
    type = nInputs;
    
    if ((mxGetM(input_pointers[0]) > (size_t)(1)) && (mxGetN(input_pointers[0]) > (size_t)(1))) {
    	std::cerr << "The original signal must be a vector!" << std::endl;
    	return;
    }
    
    if (mxGetM(input_pointers[0]) > (size_t)(1)) {
        N = mxGetM(input_pointers[0]);
    } else {
        N = mxGetN(input_pointers[0]);
    }
    
    // Normalization constant
    c = 1.0 / sqrt(N);
    
    // The original signal s and its Fourier transform F
    double *Re_s;
    double *Im_s;
    double *Re_F;
    double *Im_F;
    
    // Memory allocation
    Re_s = mxGetPr(input_pointers[0]);
    
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Re_F = mxGetPr(output_pointers[0]);
    
    output_pointers[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Im_F = mxGetPr(output_pointers[1]);
    
    // +-----------------------------+
    // | The original signal is real |
    // +-----------------------------+
    
    if (nInputs == 1) {        
        FFT_real(Re_s, Re_F, Im_F, N);
        
        for (int frequency = 0; frequency < N; ++frequency) {
            Re_F[frequency] *= c;
            Im_F[frequency] *= c;
        }
        return;
    }
    
    // +--------------------------------+
    // | The original signal is complex |
    // +--------------------------------+
    
    if (((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) || (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1])))) {
        std::cerr << "The size of the real part and the complex part are not equal!" << std::endl;
        return;
    }
    
    Im_s = mxGetPr(input_pointers[1]);
    
    FFT_complex(Re_s, Im_s, Re_F, Im_F, N);
    
    for (int frequency = 0; frequency < N; ++frequency) {
    	Re_F[frequency] *= c;
        Im_F[frequency] *= c;
    }
}