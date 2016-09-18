% Software: Fast Fourier Transform 1D and Discrete Fourier Transform 1D (for real and complex signals)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = test_dft_fft_1d()
    % Real signal
    N = 2 ^ 14;
    fprintf('1. The real signal test: N = %d\n', N);
    Signal = rand([1 N]);
    
    fprintf('Computing by the DFT 1D\n');
    [re_dft, im_dft] = dft_1d(Signal);
    
    fprintf('Computing by the FFT 1D\n');
    [re_fft, im_fft] = fft_1d(Signal);
    
    difference = sum(abs(re_dft - re_fft) + abs(im_dft - im_fft));
    fprintf('Difference: %.6f\n', difference);
    
    % Complex signal
    N = 2^12;
    fprintf('2. The complex signal test: N = %d\n', N);
    re_signal = rand([1 N]);
    im_signal = rand([1 N]);
    
    fprintf('Computing by the DFT 1D\n');
    [re_dft, im_dft] = dft_1d(re_signal, im_signal);
    
    fprintf('Computing by the FFT 1D\n');
    [re_fft, im_fft] = fft_1d(re_signal, im_signal);
    
    difference = sum(abs(re_dft - re_fft) + abs(im_dft - im_fft));
    fprintf('Difference: %.6f\n', difference);
end