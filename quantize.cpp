// Software: Quantization for JPEG Standard
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

// Block size - Constant
const int block_size = 8;

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
        std::cerr << "The number of input parameters must be exactly 2 (the matrix, and the quantization table)!" << std::endl;
        return;
    }
    
    // Number of rows
    int M = mxGetM(input_pointers[0]);
    
    // Number of columns
    int N = mxGetN(input_pointers[0]);
    
    if ((M % block_size != 0) || (N % block_size != 0)) {
        std::cerr << "The dimensions must be multiple of the block size = " << block_size << "!" << std::endl;
        return;
    }
    
    if ((mxGetM(input_pointers[1]) != block_size) || (mxGetN(input_pointers[1]) != block_size)) {
        std::cerr << "The quantization matrix dimensions are not correct!" << std::endl;
        std::cerr << "Correct sizes: " << block_size << "x" << block_size << std::endl;
        return;
    }
    
    // The matrix
    double **matrix = new double* [M];
    for (int i = 0; i < M; ++i) {
        matrix[i] = new double [N];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), M, N, matrix);
    
    // The quantization table
    double **table = new double* [block_size];
    for (int i = 0; i < block_size; ++i) {
        table[i] = new double [block_size];
    }
    
    vector2matrix(mxGetPr(input_pointers[1]), block_size, block_size, table);
    
    // Quantization
    int i = 0;
    while (i < M) {
        int j = 0;
        while (j < N) {
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    matrix[i + x][j + y] /= table[x][y];
                    matrix[i + x][j + y] = int(matrix[i + x][j + y]);
                }
            }
            j += block_size;
        }
        i += block_size;
    }
    
    // Return the quantized matrix
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(matrix, M, N, mxGetPr(output_pointers[0]));
}