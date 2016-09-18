% Software: Fast Integer Multiplication with the Fast Fourier Transform Algorithm
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_intmul_fft()
    % The number of digits of A
    nA = 2^15;
    
    % The number of digits of B
    nB = 2^15;
    
    % Randomize for number A
    A = randi([0 9], [1 nA]);
    
    % Randomize for number B
    B = randi([0 9], [1 nB]);
    
    fprintf('The number of digits of the first number: %d\n', nA);
    fprintf('The number of digits of the second number: %d\n', nB);
    
    % Fast Fourier Transform Integer Multiplication
    fprintf('1. Fast Fourier Transform Integer Multiplication\n');
    [int_fft] = intmul_fft(A, B);
    
    % Normal O(N^2) Integer Multiplication
    fprintf('2. Normal O(N^2) Integer Multiplication\n');
    [int_n2] = intmul_n2(A, B);
    
    % Check the results
    fprintf('Difference: %.6f\n', sum(abs(int_fft - int_n2)));
end