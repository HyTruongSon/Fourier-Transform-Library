// Software: Discrete Fourier Transform 1D (for real and complex signals)
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

// type = 1: The original signal is real
// type = 2: The original signal is complex
int type;

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
    if (nInputs == 0) {
        std::cerr << "Not enough input parameters!" << std::endl;
        return;
    }
    
    if (nInputs > 2) {
        std::cerr << "Maximum 2 input parameters only!" << std::endl;
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
    
    // Memory allocation
    Re_s = mxGetPr(input_pointers[0]);
    output_pointers[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    Re_F = mxGetPr(output_pointers[0]);
    
    // +-----------------------------+
    // | The original signal is real |
    // +-----------------------------+
    
    if (nInputs == 1) {
        // nOutputs = 1: Return the real part of the Fourier transform only
        // nOutputs = 2: Return both the real part and the imaginary part of the Fourier transform
        
        if (nOutputs == 1) {
            for (int k = 0; k < N; ++k) {
                Re_F[k] = 0.0;
                
                double t = 2.0 * M_PI * (double)(k) / (double)(N);
                for (int n = 0; n < N; ++n) {
                    double alpha = t * (double)(n);
                    Re_F[k] += Re_s[n] * cos(alpha);
                }
                
                Re_F[k] *= c;
            }
        } else {
            output_pointers[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
            Im_F = mxGetPr(output_pointers[1]);
        
            for (int k = 0; k < N; ++k) {
                Re_F[k] = 0.0;
                Im_F[k] = 0.0;
                
                double t = 2.0 * M_PI * (double)(k) / (double)(N);
                for (int n = 0; n < N; ++n) {
                    double alpha = t * (double)(n);
                    Re_F[k] += Re_s[n] * cos(alpha);
                    Im_F[k] -= Re_s[n] * sin(alpha);
                }
                
                Re_F[k] *= c;
                Im_F[k] *= c;
            }
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
    
    // nOutputs = 1: Return the real part of the Fourier transform only
	// nOutputs = 2: Return both the real part and the imaginary part of the Fourier transform
    
    if (nOutputs == 1) {
        for (int k = 0; k < N; ++k) {
        	Re_F[k] = 0.0;
                
        	double t = 2.0 * M_PI * (double)(k) / (double)(N);
        	for (int n = 0; n < N; ++n) {
                double alpha = t * (double)(n);
                Re_F[k] += Re_s[n] * cos(alpha) + Im_s[n] * sin(alpha);
            }
            
            Re_F[k] *= c;
        }        
    } else {
        output_pointers[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
        Im_F = mxGetPr(output_pointers[1]);
    
    	for (int k = 0; k < N; ++k) {
            Re_F[k] = 0.0;
            Im_F[k] = 0.0;
                
            double t = 2.0 * M_PI * (double)(k) / (double)(N);
        	for (int n = 0; n < N; ++n) {
                double alpha = t * (double)(n);
                double cosine = cos(alpha);
                double sine = sin(alpha);
                
                Re_F[k] += Re_s[n] * cosine + Im_s[n] * sine;
                Im_F[k] += Im_s[n] * cosine - Re_s[n] * sine;
        	}
                
            Re_F[k] *= c;
            Im_F[k] *= c;
        }
    }
}