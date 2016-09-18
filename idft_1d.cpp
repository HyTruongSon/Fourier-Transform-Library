// Software: Inverse Discrete Fourier Transform 1D (for real and complex signals)
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

// Number of time samples
int N;

// Normalization constant c = 1/sqrt(N)
double c; 

// The original signal s and its Fourier transform F
double *Re_s;
double *Im_s;
double *Re_F;
double *Im_F;

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "Exactly 2 input parameters: real and imaginary!" << std::endl;
        return;
    }
    
    if (nOutputs == 0) {
        std::cerr << "Not enough output parameters!" << std::endl;
        return;
    }
    
    if (nOutputs > 2) {
        std::cerr << "Maximum 2 output parameters only!" << std::endl;
        return;
    }
        
    if ((mxGetM(input_pointers[0]) > (size_t)(1)) && (mxGetN(input_pointers[0]) > (size_t)(1))) {
    	std::cerr << "The input signal must be a vector!" << std::endl;
    	return;
    }
    
    if ((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) && (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1]))) {
    	std::cerr << "The size of input parameters must be the same!" << std::endl;
    	return;
    }
    
    if (mxGetM(input_pointers[0]) > (size_t)(1)) {
        N = mxGetM(input_pointers[0]);
    } else {
        N = mxGetN(input_pointers[0]);
    }
    
    // Normalization constant
    c = 1.0 / sqrt(N);
    
    // Memory allocation
    Re_F = mxGetPr(input_pointers[0]);
    Im_F = mxGetPr(input_pointers[1]);
    
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Re_s = mxGetPr(output_pointers[0]);
    
    if (nOutputs == 1) {
        for (int n = 0; n < N; ++n) {
            Re_s[n] = 0;
            
            double t = 2.0 * M_PI * (double)(n) / (double)(N);
            for (int k = 0; k < N; ++k) {
                double alpha = t * (double)(k);
                Re_s[n] += Re_F[k] * cos(alpha) - Im_F[k] * sin(alpha);
            }
            
            Re_s[n] *= c;
        }
    } else {
        output_pointers[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
        Im_s = mxGetPr(output_pointers[1]);
        
        for (int n = 0; n < N; ++n) {
            Re_s[n] = 0;
            Im_s[n] = 0;
            
            double t = 2.0 * M_PI * (double)(n) / (double)(N);
            for (int k = 0; k < N; ++k) {
                double alpha = t * (double)(k);
                double cosine = cos(alpha);
                double sine = sin(alpha);
                
                Re_s[n] += Re_F[k] * cosine - Im_F[k] * sine;
                Im_s[n] += Im_F[k] * cosine + Re_F[k] * sine;
            }
            
            Re_s[n] *= c;
            Im_s[n] *= c;
        }
    }
}