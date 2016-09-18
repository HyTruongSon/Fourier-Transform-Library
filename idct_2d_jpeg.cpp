// Software: Inverse Discrete Cosine Transform 2D for JPEG
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

const int block_size = 8;
const int middle_gray = 128;

// IDCT operators
double **idct_op;

int nearestMultipleOf(int number, int block_size) {
    if (number % block_size == 0) {
        return number;
    }
    return (number / block_size + 1) * block_size;
}

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

double alpha(int u) {
    if (u == 0) {
        return 1.0 / sqrt(2.0);
    }
    return 1.0;
}

void init_IDCT_operators() {
    // Memory allocation
    int n = block_size * block_size;
    idct_op = new double* [n];
    for (int i = 0; i < n; ++i) {
        idct_op[i] = new double [n];
    }
    
    // Computation
    double c = 0.25;
    for (int x = 0; x < block_size; ++x) {
        for (int y = 0; y < block_size; ++y) {
            double t1 = (2.0 * x + 1.0) * M_PI / (2.0 * block_size);
            double t2 = (2.0 * y + 1.0) * M_PI / (2.0 * block_size);
            
            int x0 = x * block_size;
            int y0 = y * block_size;
            for (int u = 0; u < block_size; ++u) {
                for (int v = 0; v < block_size; ++v) {
                    idct_op[x0 + u][y0 + v] = c * alpha(u) * alpha(v) * cos(t1 * u) * cos(t2 * v);
                }
            }
        }
    }
}

void idct_transform(double **input_block, double **output_block) {
    for (int u = 0; u < block_size; ++u) {
        for (int v = 0; v < block_size; ++v) {
            output_block[u][v] = 0.0;
            
            int x0 = u * block_size;
            int y0 = v * block_size;
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    output_block[u][v] += input_block[x][y] * idct_op[x0 + x][y0 + y];
                }
            }
        }
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 3) {
        std::cerr << "The number of input parameters must be exactly 3 (the image, the number of rows, the number of columns)!" << std::endl;
        return;
    }
    
    // Init the operators of Inverse Discrete Cosine Transform 2D
    init_IDCT_operators();
    
    // Number of rows
    int M = mxGetScalar(input_pointers[1]);
    int nRows = nearestMultipleOf(M, block_size);
    
    // Number of columns
    int N = mxGetScalar(input_pointers[2]);
    int nCols = nearestMultipleOf(N, block_size);
    
    // DCT
    double **dct = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        dct[i] = new double [nCols];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), nRows, nCols, dct);
    
    // Inverse Discrete Cosine Transform 2D
    double **image = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        image[i] = new double [nCols];
    }
    
    double **input_block = new double* [block_size];
    double **output_block = new double* [block_size];
    
    for (int i = 0; i < block_size; ++i) {
        input_block[i] = new double [block_size];
        output_block[i] = new double [block_size];
    }
    
    int i = 0;
    while (i < nRows) {
        int j = 0;
        while (j < nCols) {
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    input_block[x][y] = dct[i + x][j + y];
                }
            }
            
            idct_transform(input_block, output_block);
            
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    image[i + x][j + y] = output_block[x][y];
                }
            }
            
            j += block_size;
        }
        i += block_size;
    }
    
    // Add back the middle of the gray scale
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            image[i][j] += middle_gray;
        }
    }
    
    // Return the IDCT
    output_pointers[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    matrix2vector(image, M, N, mxGetPr(output_pointers[0]));
}
