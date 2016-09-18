// Software: Gaussian matrix
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
        std::cerr << "The number of input parameters must be 2 (radius and variance)!" << std::endl;
        return;
    }
    
    if (nOutputs != 1) {
        std::cerr << "The number of output parameters must be 1 (the Gaussian matrix)!" << std::endl;
        return;
    }
    
    int nRows_image = mxGetM(input_pointers[0]);
    int nCols_image = mxGetN(input_pointers[0]);
    
    double **original_image = new double* [nRows_image];
    double **filtered_image = new double* [nRows_image];
    for (int row = 0; row < nRows_image; ++row) {
        original_image[row] = new double [nCols_image];
        filtered_image[row] = new double [nCols_image];
    }
    
    vector2matrix(mxGetPr(input_pointers[0]), nRows_image, nCols_image, original_image);
    
    int nRows_filter = mxGetM(input_pointers[1]);
    int nCols_filter = mxGetN(input_pointers[1]);
    
    double **filter = new double* [nRows_filter];
    for (int row = 0; row < nRows_filter; ++row) {
        filter[row] = new double [nCols_filter];
    }
    
    vector2matrix(mxGetPr(input_pointers[1]), nRows_filter, nCols_filter, filter);
    
    int x0 = nRows_filter / 2;
    int y0 = nCols_filter / 2;
    
    int x1 = nRows_image - nRows_filter + x0;
    int y1 = nCols_image - nCols_filter + y0;
    
    for (int i = 0; i < nRows_image; ++i) {
        for (int j = 0; j < nCols_image; ++j) {
            if ((i >= x0) && (i <= x1) && (j >= y0) && (j <= y1)) {
                int x = i - x0;
                int y = j - y0;
                double sum = 0.0;
                for (int u = 0; u < nRows_filter; ++u) {
                    for (int v = 0; v < nCols_filter; ++v) {
                        sum += original_image[x + u][y + v] * filter[u][v];
                    }
                }
                filtered_image[i][j] = (int)(sum);
            } else {
                filtered_image[i][j] = original_image[i][j];
            }
        }
    }
    
    output_pointers[0] = mxCreateDoubleMatrix(nRows_image, nCols_image, mxREAL);
    matrix2vector(filtered_image, nRows_image, nCols_image, mxGetPr(output_pointers[0]));
}