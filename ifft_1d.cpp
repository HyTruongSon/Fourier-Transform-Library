// Software: Inverse Fast Fourier Transform 1D (for real and complex signals)
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time complexity: O(NlogN)
// Space complexity: O(N).

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

// Fast Fourier Transform 1D (with the different harmonic coefficient) for complex signals
// It becomes the Inverse Fast Fourier Transform 1D
void iFFT(double *Re_Signal, double *Im_Signal, double *Re_F, double *Im_F, int N, int t) {
        if (N == 1) {
        Re_F[0] = Re_Signal[0];
        Im_F[0] = Im_Signal[0];
        return;
    }
    
    int half = N / 2;
    iFFT(Re_Signal, Im_Signal, Re_F, Im_F, half, 2 * t);
    iFFT(Re_Signal + t, Im_Signal + t, Re_F + half, Im_F + half, half, 2 * t);
    
    for (int k = 0; k < half; ++k) {
        double r1 = Re_F[k];
        double i1 = Im_F[k];
        
        double r2 = Re_F[k + half];
        double i2 = Im_F[k + half];
        
        double alpha = 2.0 * M_PI * (double)(k) / (double)(N);
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
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be 2 (real and imaginary parts)!" << std::endl;
        return;
    }
    
    if (nOutputs == 0) {
        std::cerr << "Not enough output parameters!" << std::endl;
        return;
    }
    
    if (nOutputs > 2) {
        std::cerr << "Maximum 2 output parameters!" << std::endl;
        return;
    }
    
    type = nOutputs;
    
    if ((mxGetM(input_pointers[0]) > (size_t)(1)) && (mxGetN(input_pointers[0]) > (size_t)(1))) {
    	std::cerr << "The original signal must be a vector!" << std::endl;
    	return;
    }
    
    if ((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) || (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1]))) {
    	std::cerr << "The size of two parameters must be the same!" << std::endl;
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
    double *Re_Signal;
    double *Im_Signal;
    double *Re_F;
    double *Im_F;
    
    // Memory allocation
    Re_F = mxGetPr(input_pointers[0]);
    Im_F = mxGetPr(input_pointers[1]);
    
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Re_Signal = mxGetPr(output_pointers[0]);
    
    if (nOutputs == 2) {
        output_pointers[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
        Im_Signal = mxGetPr(output_pointers[1]);
    } else {
        Im_Signal = new double [N];
    }
    
    // The order of parameters change compared to FFT
    iFFT(Re_F, Im_F, Re_Signal, Im_Signal, N, 1);
    
    for (int frequency = 0; frequency < N; ++frequency) {
    	Re_Signal[frequency] *= c;
        Im_Signal[frequency] *= c;
    }
}