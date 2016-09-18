// Software: Fast Polynomial Multiplication with the Fast Fourier Transform Algorithm
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

int nearest_power_2(int m) {
    int k = int(log(m) / log(2));
    int n = (1 << k);
    if (n >= m) {
        return n;
    }
    return 2 * n;
}

void FFT(double *P, double *Re_F, double *Im_F, int N, int t, double Re_w, double Im_w) {
    if (N == 1) {
        Re_F[0] = P[0];
        Im_F[0] = 0.0;
        return;
    }
    
    double Re_w2 = Re_w * Re_w - Im_w * Im_w;
    double Im_w2 = 2 * Re_w * Im_w;
    
    int half = N / 2;
    FFT(P, Re_F, Im_F, half, 2 * t, Re_w2, Im_w2);
    FFT(P + t, Re_F + half, Im_F + half, half, 2 * t, Re_w2, Im_w2);
    
    double Re_x = 1.0;
    double Im_x = 0.0;
    
    for (int k = 0; k < half; ++k) {
        double Re_E = Re_F[k];
        double Im_E = Im_F[k];
        
        double Re_O = Re_F[k + half];
        double Im_O = Im_F[k + half];
        
        double r = Re_x * Re_O - Im_x * Im_O;
        double i = Im_x * Re_O + Re_x * Im_O;
        
        Re_F[k] = Re_E + r;
        Im_F[k] = Im_E + i;
        
        Re_F[k + half] = Re_E - r;
        Im_F[k + half] = Im_E - i;
        
        double Re_next = Re_x * Re_w - Im_x * Im_w;
        double Im_next = Im_x * Re_w + Re_x * Im_w;
        
        Re_x = Re_next;
        Im_x = Im_next;
    }
}

void iFFT(double *Re_P, double *Im_P, double *Re_F, double *Im_F, int N, int t, double Re_w, double Im_w) {
    if (N == 1) {
        Re_F[0] = Re_P[0];
        Im_F[0] = Im_P[0];
        return;
    }
    
    double Re_w2 = Re_w * Re_w - Im_w * Im_w;
    double Im_w2 = 2 * Re_w * Im_w;
    
    int half = N / 2;
    iFFT(Re_P, Im_P, Re_F, Im_F, half, 2 * t, Re_w2, Im_w2);
    iFFT(Re_P + t, Im_P + t, Re_F + half, Im_F + half, half, 2 * t, Re_w2, Im_w2);
    
    double Re_x = 1.0;
    double Im_x = 0.0;
    
    for (int k = 0; k < half; ++k) {
        double Re_E = Re_F[k];
        double Im_E = Im_F[k];
        
        double Re_O = Re_F[k + half];
        double Im_O = Im_F[k + half];
        
        double r = Re_x * Re_O - Im_x * Im_O;
        double i = Im_x * Re_O + Re_x * Im_O;
        
        Re_F[k] = Re_E + r;
        Im_F[k] = Im_E + i;
        
        Re_F[k + half] = Re_E - r;
        Im_F[k + half] = Im_E - i;
        
        double Re_next = Re_x * Re_w - Im_x * Im_w;
        double Im_next = Im_x * Re_w + Re_x * Im_w;
        
        Re_x = Re_next;
        Im_x = Im_next;
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2 (polynomials)!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) > size_t(1)) && (mxGetN(input_pointers[0]) > size_t(1))) {
        std::cerr << "The input parameters must be vectors!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[1]) > size_t(1)) && (mxGetN(input_pointers[1]) > size_t(1))) {
        std::cerr << "The input parameters must be vectors!" << std::endl;
        return;
    }
    
    int nA, nB, nC, N;
    double *A, *B, *C, *Re_C, *Im_C, *temp;
    double *Re_FA, *Im_FA, *Re_FB, *Im_FB, *Re_FC, *Im_FC;
    
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
    N = nearest_power_2(nC);
    
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
    
    temp = mxGetPr(input_pointers[0]);
    for (int i = 0; i < nA; ++i) {
        A[i] = temp[i];
    }
    for (int i = nA; i < N; ++i) {
        A[i] = 0.0;
    }
    
    temp = mxGetPr(input_pointers[1]);
    for (int i = 0; i < nB; ++i) {
        B[i] = temp[i];
    }
    for (int i = nB; i < N; ++i) {
        B[i] = 0.0;
    }
    
    // Fast Fourier Transform 1D
    double alpha = 2.0 * M_PI / double(N);
    double Re_w1 = cos(alpha);
    double Im_w1 = sin(alpha);
    
    FFT(A, Re_FA, Im_FA, N, 1, Re_w1, Im_w1);
    FFT(B, Re_FB, Im_FB, N, 1, Re_w1, Im_w1);
    
    for (int i = 0; i < N; ++i) {
        Re_FC[i] = Re_FA[i] * Re_FB[i] - Im_FA[i] * Im_FB[i];
        Im_FC[i] = Re_FA[i] * Im_FB[i] + Im_FA[i] * Re_FB[i];
    }
    
    // Inverse Fast Fourier Transform 1D
    double Re_w2 = cos(alpha);
    double Im_w2 = - sin(alpha);
    
    iFFT(Re_FC, Im_FC, Re_C, Im_C, N, 1, Re_w2, Im_w2);
    
    // Return the multiplication polynomial
    output_pointers[0] = mxCreateDoubleMatrix(nC, 1, mxREAL);
    C = mxGetPr(output_pointers[0]);
    
    for (int i = 0; i < nC; ++i) {
        C[i] = Re_C[i] / double(N);
    }
}