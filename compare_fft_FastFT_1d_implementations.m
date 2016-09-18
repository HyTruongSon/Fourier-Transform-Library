% Software: Fast Fourier Transform 1D
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = compare_fft_FastFT_1d_implementations()
    fprintf('Comparison the performance of fft and FastFT 1D:');
    N = 2 ^ 20;
    fprintf('- Generating the signals: N = %d\n', N);
    signal = rand([1 N]);
    
    fprintf('- fft_1d and ifft_1d\n');
    fprintf('  Time complexity: O(NlogN)\n');
    fprintf('  Space complexity: O(N)\n');
    
    [re1, im1] = fft_1d(signal);
    [signal1] = ifft_1d(re1, im1);
    
    fprintf('  Difference: %.6f\n', sum(abs(signal - signal1')));
    
    fprintf('- FastFT_1d and iFastFT_1d\n');
    fprintf('  Time complexity: O(NlogN)\n');
    fprintf('  Space complexity: O(NlogN)\n');
    
    [re2, im2] = FastFT_1d(signal);
    [signal2] = iFastFT_1d(re2, im2);
    
    fprintf('  Difference: %.6f\n', sum(abs(signal - signal2')));    
end