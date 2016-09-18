% Software: Fast Polynomial Multiplication with the Fast Fourier Transform Algorithm
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_polymul_fft()
    % Constants
    nA = 2^16; % First polynomial degree
    nB = 2^16; % Second polynomial degree
    
    % Input polynomial randomization
    A = rand([1 nA]);
    B = rand([1 nB]);
    
    fprintf('First polynomial degree: %d\n', nA);
    fprintf('Second polynomial degree: %d\n', nB);
    
    % Fast Fourier Transform Polynomial Multiplication
    fprintf('1. Fast Fourier Transform Polynomial Multiplication\n');
    [poly_fft] = polymul_fft(A, B);
    
    % Normal O(N^2) Polynomial Multiplication
    fprintf('2. Normal O(N^2) Polynomial Multiplication\n');
    [poly_n2] = polymul_n2(A, B);
    
    % Check the results
    fprintf('Difference: %.6f\n', sum(abs(poly_fft - poly_n2)));
end