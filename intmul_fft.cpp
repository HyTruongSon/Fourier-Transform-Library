// Software: Fast Integer Multiplication (Fast Fourier Transform)
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time Complexity: O(NlogN)
// Space Complexity: O(N)

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

int nearestPower2(int n) {
    int N = (1 << (int)(log(n) / log(2)));
    if (N >= n) {
        return N;
    }
    return 2 * N;
}

int nearestInteger(double num) {
    double i = (int)(num);
    
    if (num == i) {
        return i;
    }
    
    if (num > i) {
        if (num - i > 0.5) {
            return i + 1;
        }
        return i;
    }
    
    if (i - num > 0.5) {
        return i - 1;
    }
    return i;
}

// Fast Fourier Transform (for real signals)
void FFT_real(double *P, double *Re_FP, double *Im_FP, int N, int t, double Re_w, double Im_w) {
    if (N == 1) {
        Re_FP[0] = P[0];
        Im_FP[0] = 0.0;
        return;
    }
    
    double Re_w2 = Re_w * Re_w - Im_w * Im_w;
    double Im_w2 = 2.0 * Re_w * Im_w;
    
    int half = N / 2;
    FFT_real(P, Re_FP, Im_FP, half, 2 * t, Re_w2, Im_w2);
    FFT_real(P + t, Re_FP + half, Im_FP + half, half, 2 * t, Re_w2, Im_w2);
    
    double Re_x = 1.0;
    double Im_x = 0.0;
    
    for (int j = 0; j < half; ++j) {
        double Re_E = Re_FP[j];
        double Im_E = Im_FP[j];
        
        double Re_O = Re_FP[j + half];
        double Im_O = Im_FP[j + half];
        
        double r = Re_x * Re_O - Im_x * Im_O;
        double i = Re_x * Im_O + Im_x * Re_O;
        
        Re_FP[j] = Re_E + r;
        Im_FP[j] = Im_E + i;
        
        Re_FP[j + half] = Re_E - r;
        Im_FP[j + half] = Im_E - i;
        
        double Re_nextx = Re_x * Re_w - Im_x * Im_w;
        double Im_nextx = Re_x * Im_w + Im_x * Re_w;
        
        Re_x = Re_nextx;
        Im_x = Im_nextx;
    }
}

// Fast Fourier Transform (for complex signals)
void FFT_complex(double *Re_P, double *Im_P, double *Re_FP, double *Im_FP, int N, int t, double Re_w, double Im_w) {
    if (N == 1) {
        Re_FP[0] = Re_P[0];
        Im_FP[0] = Im_P[0];
        return;
    }
    
    double Re_w2 = Re_w * Re_w - Im_w * Im_w;
    double Im_w2 = 2.0 * Re_w * Im_w;
    
    int half = N / 2;
    FFT_complex(Re_P, Im_P, Re_FP, Im_FP, half, 2 * t, Re_w2, Im_w2);
    FFT_complex(Re_P + t, Im_P + t, Re_FP + half, Im_FP + half, half, 2 * t, Re_w2, Im_w2);
    
    double Re_x = 1.0;
    double Im_x = 0.0;
    
    for (int j = 0; j < half; ++j) {
        double Re_E = Re_FP[j];
        double Im_E = Im_FP[j];
        
        double Re_O = Re_FP[j + half];
        double Im_O = Im_FP[j + half];
        
        double r = Re_x * Re_O - Im_x * Im_O;
        double i = Re_x * Im_O + Im_x * Re_O;
        
        Re_FP[j] = Re_E + r;
        Im_FP[j] = Im_E + i;
        
        Re_FP[j + half] = Re_E - r;
        Im_FP[j + half] = Im_E - i;
        
        double Re_nextx = Re_x * Re_w - Im_x * Im_w;
        double Im_nextx = Re_x * Im_w + Im_x * Re_w;
        
        Re_x = Re_nextx;
        Im_x = Im_nextx;
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) > size_t(1)) && (mxGetN(input_pointers[0]) > size_t(1))) {
        std::cerr << "The first input parameter must be a vector!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[1]) > size_t(1)) && (mxGetN(input_pointers[1]) > size_t(1))) {
        std::cerr << "The second input parameter must be a vector!" << std::endl;
        return;
    }
    
    int nA, nB, nC, N;
    double *A, *B, *Re_C, *Im_C;
    double *Re_FA, *Im_FA, *Re_FB, *Im_FB, *Re_FC, *Im_FC;
    int *C;
    
    if (mxGetM(input_pointers[0]) > size_t(1)) {
        nA = mxGetM(input_pointers[0]);
    } else {
        nA = mxGetN(input_pointers[0]);
    }
    
    if (mxGetM(input_pointers[1]) > size_t(1)) {
        nB = mxGetM(input_pointers[1]);
    } else {
        nB = mxGetN(input_pointers[1]);
    }
    
    nC = nA + nB - 1;
    N = nearestPower2(nC);
    
    // Memory allocation
    A = new double [N];
    B = new double [N];
    Re_C = new double [N];
    Im_C = new double [N];
    Re_FA = new double [N];
    Im_FA = new double [N];
    Re_FB = new double [N];
    Im_FB = new double [N];
    Re_FC = new double [N];
    Im_FC = new double [N];
    C = new int [N];
    
    // Get the first integer (polynomial coefficients)
    double *inA = mxGetPr(input_pointers[0]);
    for (int i = 0; i < nA; ++i) {
        A[i] = inA[nA - i - 1];
    }
    for (int i = nA; i < N; ++i) {
        A[i] = 0.0;
    }
    
    // Get the second integer (polynomial coefficients)
    double *inB = mxGetPr(input_pointers[1]);
    for (int i = 0; i < nB; ++i) {
        B[i] = inB[nB - i - 1];
    }
    for (int i = nB; i < N; ++i) {
        B[i] = 0.0;
    }
    
    // Fast Fourier Transform
    double alpha = 2.0 * M_PI / (double)(N);
    double Re_w = cos(alpha);
    double Im_w = sin(alpha);
    
    FFT_real(A, Re_FA, Im_FA, N, 1, Re_w, Im_w);
    FFT_real(B, Re_FB, Im_FB, N, 1, Re_w, Im_w);
    
    for (int i = 0; i < N; ++i) {
        Re_FC[i] = Re_FA[i] * Re_FB[i] - Im_FA[i] * Im_FB[i];
        Im_FC[i] = Re_FA[i] * Im_FB[i] + Im_FA[i] * Re_FB[i];
    }
    
    // Inverse Fast Fourier Transform
    FFT_complex(Re_FC, Im_FC, Re_C, Im_C, N, 1, Re_w, -Im_w);
    
    for (int i = 0; i < nC; ++i) {
        C[i] = nearestInteger(Re_C[i] / (double)(N));
    }
    
    // Hornell scheme
    int *P = new int [nA + nB + 1];
    int *T = new int [nA + nB + 1];
    int nP, nT;
    
    nP = 0;
    int c = C[nC - 1];
    while (c > 0) {
        P[nP] = c % 10;
        c /= 10;
        ++nP;
    }
    
    for (int i = nC - 1; i > 0; --i) {
        // Multiplies with 10
        for (int j = nP - 1; j >= 0; --j) {
            P[j + 1] = P[j];
        }
        ++nP;
        P[0] = 0;
        
        // Get a new coefficient
        nT = 0;
        c = C[i - 1];
        while (c > 0) {
            T[nT] = c % 10;
            c /= 10;
            ++nT;
        }
        
        // Sum up
        int n = max(nP, nT);
        int r = 0;
        
        for (int j = 0; j < n; ++j) {
            int u = 0;
            if (j < nP) {
                u = P[j];
            }
            
            int v = 0;
            if (j < nT) {
                v = T[j];
            }
            
            P[j] = u + v + r;
            r = P[j] / 10;
            P[j] %= 10;
        }
        
        if (r > 0) {
            nP = n + 1;
            P[nP - 1] = r;
        } else {
            nP = n;
        }
    }
    
    // Return the result
    output_pointers[0] = mxCreateDoubleMatrix(nP, 1, mxREAL);
    double *result = mxGetPr(output_pointers[0]);
    
    for (int i = 0; i < nP; ++i) {
        result[i] = P[nP - i - 1];
    }
}