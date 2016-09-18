// Software: Discrete Cosine Transform 2D for JPEG
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

// DCT operators
double **dct_op;

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

void init_DCT_operators() {
    // Memory allocation
    int n = block_size * block_size;
    dct_op = new double* [n];
    for (int i = 0; i < n; ++i) {
        dct_op[i] = new double [n];
    }
    
    // Computation
    for (int u = 0; u < block_size; ++u) {
        for (int v = 0; v < block_size; ++v) {
            double c = 0.25 * alpha(u) * alpha(v);
            double t1 = (double)(u) * M_PI / (2.0 * (double)(block_size));
            double t2 = (double)(v) * M_PI / (2.0 * (double)(block_size));
            
            int x0 = u * block_size;
            int y0 = v * block_size;
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    dct_op[x0 + x][y0 + y] = c * cos((2.0 * x + 1.0) * t1) * cos((2.0 * y + 1.0) * t2);
                }
            }
        }
    }
}

void dct_transform(double **input_block, double **output_block) {
    for (int u = 0; u < block_size; ++u) {
        for (int v = 0; v < block_size; ++v) {
            output_block[u][v] = 0.0;
            
            int x0 = u * block_size;
            int y0 = v * block_size;
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    output_block[u][v] += input_block[x][y] * dct_op[x0 + x][y0 + y];
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
    
    // Init the operators of Discrete Cosine Transform 2D
    init_DCT_operators();
    
    // Number of rows
    int M = mxGetScalar(input_pointers[1]);
    int nRows = nearestMultipleOf(M, block_size);
    
    // Number of columns
    int N = mxGetScalar(input_pointers[2]);
    int nCols = nearestMultipleOf(N, block_size);
    
    // Image
    double **image = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        image[i] = new double [nCols];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), M, N, image);
    
    // Add more rows and columns such that the sizes must be multiple of the block size
    for (int i = M; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            image[i][j] = 0.0;
        }
    }
    
    for (int j = N; j < nCols; ++j) {
        for (int i = 0; i < nRows; ++i) {
            image[i][j] = 0.0;
        }
    }
    
    // Subtract to the middle of the gray scale [0..255]
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            image[i][j] -= middle_gray;
        }
    }
    
    // Discrete Cosine Transform 2D
    double **dct = new double* [nRows];
    for (int i = 0; i < nRows; ++i) {
        dct[i] = new double [nCols];
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
                    input_block[x][y] = image[i + x][j + y];
                }
            }
            
            dct_transform(input_block, output_block);
            
            for (int x = 0; x < block_size; ++x) {
                for (int y = 0; y < block_size; ++y) {
                    dct[i + x][j + y] = output_block[x][y];
                }
            }
            
            j += block_size;
        }
        i += block_size;
    }
    
    // Return the DCT
    output_pointers[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    matrix2vector(dct, nRows, nCols, mxGetPr(output_pointers[0]));
}
