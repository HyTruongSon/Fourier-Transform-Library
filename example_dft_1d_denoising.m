% Software: Discrete Fourier Transform 1D (for real and complex signals)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_dft_1d_denoising()
    % Data generation
    t = [0 : 0.1 : 127.9];
    N = size(t, 2);
    original_signal = 5 * sin(t / 10) + 2 * cos(t / 4) + sin(2 * t + pi / 2);
    noise = rand([1 N]) - 0.5;
    noisy_signal = original_signal + noise;
    
    figure(1);
    plot(t(:), noisy_signal(:), 'r-');
    hold on;
    plot(t(:), original_signal(:), 'b-');
    title('The original signal and the noisy signal');
    xlabel('Time');
    ylabel('Magnitude');
    
    % DFT 1D
    [original_re, original_im] = dft_1d(original_signal);
    [original_amp] = amplitude_1d(original_re, original_im);
    
    [noisy_re, noisy_im] = dft_1d(noisy_signal);
    [noisy_amp] = amplitude_1d(noisy_re, noisy_im);
    
    figure(2);
    plot(t(:), original_amp(:), 'b.');
    hold on;
    title('The amplitude of the original signal in the frequency domain');
    xlabel('Time');
    ylabel('Amplitude');
    
    figure(3);
    plot(t(:), noisy_amp(:), 'r.');
    hold on;
    title('The amplitude of the noisy signal in the frequency domain');
    xlabel('Time');
    ylabel('Amplitude');
    
    % Denoising
    % threshold = sum(noisy_amp) / N;
    threshold = 1.0;
    
    denoised_re = noisy_re;
    denoised_im = noisy_im;
    denoised_re(noisy_amp < threshold) = 0.0;
    denoised_im(noisy_amp < threshold) = 0.0;
    denoised_amp = amplitude_1d(denoised_re, denoised_im);
    
    figure(4);
    plot(t(:), denoised_amp(:), 'g.');
    hold on;
    title('The amplitude of the denoised signal in the frequency domain');
    xlabel('Time');
    ylabel('Amplitude');
    
    % Recovery
    [denoised_signal] = idft_1d(denoised_re, denoised_im);
    
    figure(5);
    plot(t(:), denoised_signal(:), 'g-');
    hold on;
    plot(t(:), original_signal(:), 'b-');
    title('The original signal and the noisy signal');
    xlabel('Time');
    ylabel('Magnitude');
end