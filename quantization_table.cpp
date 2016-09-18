// Software: Quantization Tables for JPEG Standard
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

// Quantization table - Luminance for JPEG Standard - ID = 1
const int luminance[block_size][block_size] = {
    {16,    11,     10,     16,     24,     40,     51,     61},
    {12,    12,     14,     19,     26,     58,     60,     55},
    {14,    13,     16,     24,     40,     57,     69,     56},
    {14,    17,     22,     29,     51,     87,     80,     62},
    {18,    22,     37,     56,     68,     109,    103,    77},
    {24,    35,     55,     64,     81,     104,    113,    92},
    {49,    64,     78,     87,     103,    121,    120,    101},
    {72,    92,     95,     98,     112,    100,    103,    99}
};

// Quantization table - Chrominance for JPEG Standard - ID = 2
const int chrominance[block_size][block_size] = {
    {17,    18,     24,     47,     99,     99,     99,     99},
    {18,    21,     26,     66,     99,     99,     99,     99},
    {24,    26,     56,     99,     99,     99,     99,     99},
    {47,    66,     99,     99,     99,     99,     99,     99},
    {99,    99,     99,     99,     99,     99,     99,     99},
    {99,    99,     99,     99,     99,     99,     99,     99},
    {99,    99,     99,     99,     99,     99,     99,     99},
    {99,    99,     99,     99,     99,     99,     99,     99}
};

void matrix2vector(const int input[block_size][block_size], double *output) {
    for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
            output[j * block_size + i] = input[i][j];
        }
    }
}

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 1) {
        std::cerr << "The number of parameters must be exactly 1 (the id of the table)!" << std::endl;
        return;
    }
    
    int id = mxGetScalar(input_pointers[0]);
    
    if ((id < 1) || (id > 2)) {
        std::cerr << "Unknown id!" << std::endl;
        std::cerr << "ID = 1: Luminance" << std::endl;
        std::cerr << "ID = 2: Chrominance" << std::endl; 
        return;
    }
    
    // Return the quantization table
    output_pointers[0] = mxCreateDoubleMatrix(block_size, block_size, mxREAL);
    
    // Luminance quantization table
    if (id == 1) {
        matrix2vector(luminance, mxGetPr(output_pointers[0]));
    }
    
    // Chrominance quantization table
    if (id == 2) {
        matrix2vector(chrominance, mxGetPr(output_pointers[0]));
    }
}