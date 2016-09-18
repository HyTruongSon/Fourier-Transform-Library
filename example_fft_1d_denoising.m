% Software: Fast Fourier Transform 1D (for real and complex signals)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_1d_denoising()
    % Generate a random signal
    t = [0 : 0.1 : 20 * pi];
    t = t(1:2^9);
    c = 5.0;
    sine1 = c * sin(t / 10);
    sine2 = c * sin(t / 2);
    sine3 = c * sin(t / 3);
    sine4 = c * sin(t / 6 + pi / 4);
    cosine1 = c * cos(t / 20);
    cosine2 = c * cos(t / 2);
    cosine3 = c * cos(t / 3 + pi / 2.5);
    cosine4 = c * cos(t / 5 - pi / 3);
    original_signal = sine1 + sine2 + sine3 + sine4 + cosine1 + cosine2 + cosine3 + cosine4;
    N = size(t, 2);
    noise = 2.0 * rand([1 N]) - 1.0;
    noisy_signal = original_signal + noise;
    
    figure(1);
    plot(t(:), noisy_signal(:), 'r-');
    hold on;
    plot(t(:), original_signal(:), 'b-');
    title('The original signal and the noisy signal');
    xlabel('Time');
    ylabel('Magnitude');
    
    % Fast Fourier Transform 1D
    [original_re, original_im] = fft_1d(original_signal);
    [original_amp] = amplitude_1d(original_re, original_im);
    
    [noisy_re, noisy_im] = fft_1d(noisy_signal);
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
    threshold = 4.0;
    
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
    
    % Recovery by Fast Fourier Tranform
    [denoised_signal] = ifft_1d(denoised_re, denoised_im);
    
    figure(5);
    plot(t(:), denoised_signal(:), 'g-');
    hold on;
    plot(t(:), original_signal(:), 'b-');
    title('The original signal and the noisy signal');
    xlabel('Time');
    ylabel('Magnitude');
end