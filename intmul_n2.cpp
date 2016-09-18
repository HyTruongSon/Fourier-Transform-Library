// Software: Normal Integer Multiplication
// Author: Hy Truong Son
// Position: PhD Student
// Institution: Department of Computer Science, The University of Chicago
// Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
// Website: http://people.inf.elte.hu/hytruongson/
// Copyright 2016 (c) Hy Truong Son. All rights reserved.

// Time Complexity: O(N^2)
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

void mexFunction(int nOutputs, mxArray *output_pointers[], int nInputs, const mxArray *input_pointers[]) {
    if (nInputs != 2) {
        std::cerr << "The number of input parameters must be exactly 2!" << std::endl;
        return;
    }
    
    if (nOutputs != 1) {
        std::cerr << "The number of output parameters must be exactly 1!" << std::endl;
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
    
    int N, M;
    
    if (mxGetM(input_pointers[0]) > size_t(1)) {
        M = mxGetM(input_pointers[0]);
    } else {
        M = mxGetN(input_pointers[0]);
    }
    
    if (mxGetM(input_pointers[1]) > size_t(1)) {
        N = mxGetM(input_pointers[1]);
    } else {
        N = mxGetN(input_pointers[1]);
    }
    
    double *in1 = mxGetPr(input_pointers[0]);
    double *in2 = mxGetPr(input_pointers[1]);
    
    int *a = new int [M];    
    for (int i = 0; i < M; ++i) {
        a[M - i - 1] = int(in1[i]);
    }
    
	int *b = new int [N];
    for (int i = 0; i < N; ++i) {
        b[N - i - 1] = int(in2[i]);
    }
    
    int *c = new int [N + M + 1];
    int *temp = new int [N + M + 1];
    
    int K = 0;
    int T;
    
    for (int i = 0; i < N; ++i) {
        int r = 0;
        for (int j = 0; j < M; ++j) {
            temp[j] = b[i] * a[j] + r;
            r = temp[j] / 10;
            temp[j] %= 10;
        }
        
        if (r > 0) {
            T = M + 1;
            temp[T - 1] = r;
        } else {
            T = M;
        }
        
        if (i == 0) {
            K = T;
            for (int j = 0; j < K; ++j) {
                c[j] = temp[j];
            }
        } else {
            for (int j = T - 1; j >= 0; --j) {
                temp[j + i] = temp[j];
            }
            for (int j = 0; j < i; ++j) {
                temp[j] = 0;
            }
            T += i;
            
            r = 0;
            for (int j = 0; j < max(T, K); ++j) {
                int u = 0;
                if (j < T) {
                    u = temp[j];
                }
                int v = 0;
                if (j < K) {
                    v = c[j];
                }
                c[j] = u + v + r;
                r = c[j] / 10;
                c[j] %= 10;
            }
            if (r > 0) {
                K = max(T, K) + 1;
                c[K - 1] = r;
            } else {
                K = max(T, K);
            }
        }
    }
    
    output_pointers[0] = mxCreateDoubleMatrix(K, 1, mxREAL);
    double *result = mxGetPr(output_pointers[0]);
    for (int i = 0; i < K; ++i) {
        result[i] = c[K - i - 1];
    }
}