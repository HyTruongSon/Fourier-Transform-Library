% Software: Fast Fourier Transform 1D (for real and complex signals)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_1d()
    % Generate a random signal
    t = [0 : 0.1 : 200 * pi];
    t = t(1:2^12);
    c = 5.0;
    sine1 = c * sin(t / 10);
    sine2 = c * sin(t / 2);
    sine3 = c * sin(t / 3);
    sine4 = c * sin(t / 6 + pi / 4);
    cosine1 = c * cos(t / 20);
    cosine2 = c * cos(t / 2);
    cosine3 = c * cos(t / 3 + pi / 2.5);
    cosine4 = c * cos(t / 5 - pi / 3);
    signal = sine1 + sine2 + sine3 + sine4 + cosine1 + cosine2 + cosine3 + cosine4;
    
    % Plot the original signal
    figure(1);
    title('The original signal');
    xlabel('Time');
    ylabel('Magnitude');
    hold on;
    plot(t(:), signal(:), 'b-');
    
    % Fast Fourier Transform
    [re, im] = fft_1d(signal);
    amp = amplitude_1d(re, im);
    phase = phase_1d(re, im);
    
    % Plot the amplitude in the frequency domain
    figure(2);
    title('The amplitude in the frequency domain after DFT');
    xlabel('Time');
    ylabel('Amplitude');
    hold on;
    plot(1:size(amp, 1), amp(:), 'g.');
    
    % Inverse FFT 1D
    [r, i] = ifft_1d(re, im);
    
    % Plot the reversed signal (real part) and the imaginary part (= 0)
    figure(3);
    title('The reversed signal from Inverse DFT');
    xlabel('Time');
    ylabel('Magnitude');
    hold on;
    plot(t(:), r(:), 'r-');
    hold on;
    plot(t(:), i(:), 'g-');
end