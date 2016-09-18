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

int M, N;
double **Re;
double **Im;
double **Phase;

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
        std::cerr << "Exactly 2 input parameters!" << std::endl;
        return;
    }
    
    if (nOutputs != 1) {
        std::cerr << "Exactly 1 output parameters!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[0]) != mxGetM(input_pointers[1])) || (mxGetN(input_pointers[0]) != mxGetN(input_pointers[1]))) {
        std::cerr << "The size of the real part and the imaginary part must be the same!" << std::endl;
        return;
    }
    
    int M = mxGetM(input_pointers[0]);
    int N = mxGetN(input_pointers[0]);
    
    Re = new double* [M];
    Im = new double* [N];
    Phase = new double* [M];
    
    for (int row = 0; row < M; ++row) {
        Re[row] = new double [N];
        Im[row] = new double [N];
        Phase[row] = new double [N];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), M, N, Re);
    vector2matrix(mxGetPr(input_pointers[1]), M, N, Im);
    
    for (int row = 0; row < M; ++row) {
        for (int column = 0; column < N; ++column) {
            Phase[row][column] = atan2(Im[row][column], Re[row][column]);
        }
    }
    
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(Phase, M, N, mxGetPr(output_pointers[0]));
}