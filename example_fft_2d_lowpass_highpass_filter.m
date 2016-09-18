% Software: Fast Fourier Transform 2D for Image Processing
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_2d_lowpass_highpass_filter()
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    
    % Constants
    radius = 50; % Pixels
    
    original_img = imread(img_fn);
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    if (size(original_img, 3) == 3)
        original_img = rgb2gray(imread(img_fn));
    end
    
    % Fast Fourier Transform 2D of the original image
    [re, im] = fft_2d(double(original_img));
    
    % Compute amplitude of the original image
    [visual_amplitude] = visual_amplitude_2d(re, im); % Only for visualization
    
    figure(2);
    title('The visual amplitude in the frequency domain (original image)');
    hold on;
    imshow(uint8(visual_amplitude));
    
    % Compute the circle for Lowpass filter
    N = size(original_img, 1);
    active = zeros(N, N);
    x0 = N / 2;
    y0 = N / 2;
    for x = x0 - radius : x0 + radius
        delta = radius * radius - (x - x0) * (x - x0);
        if delta >= 0
            r = int32(sqrt(delta));
            active(x, y0 - r : y0 + r) = 1;
        end
    end
    
    % Zeroing out the outside of the circle
    lowpass_amplitude = visual_amplitude;
    highpass_amplitude = visual_amplitude;
    
    lowpass_amplitude(active == 0) = 0.0;
    highpass_amplitude(active == 1) = 0.0;
    
    % Note that: the visual amplitude is shifted comparing to the original
    % (correct) amplitude which is corresponded to re and im matrices.
    
    active = zeros(N, N);
    for x = 1 : min(radius, N)
        delta = int32(sqrt(radius * radius - x * x));
        active(x, 1 : 1 + delta) = 1;
        active(x, N - delta : N) = 1;
    end
    
    for x = max(1, N - radius) : N
        delta = int32(sqrt(radius * radius - (x - N) * (x - N)));
        active(x, 1 : 1 + delta) = 1;
        active(x, N - delta : N) = 1;
    end
    
    % Lowpass filter
    lowpass_re = re;
    lowpass_im = im;
    
    lowpass_re(active == 0) = 0.0;
    lowpass_im(active == 0) = 0.0;
    
    % Highpass filter
    highpass_re = re;
    highpass_im = im;
    
    highpass_re(active == 1) = 0.0;
    highpass_im(active == 1) = 0.0;
    
    % Reconstruct the image
    [lowpass_img] = ifft_2d(lowpass_re, lowpass_im);
    [highpass_img] = ifft_2d(highpass_re, highpass_im);
    
    figure(3);
    title('The visual amplitude in the lowpass-filtered frequency domain');
    hold on;
    imshow(uint8(lowpass_amplitude));
    
    figure(4);
    title('The lowpass-filter in the original Fourier frequency space');
    hold on;
    imshow(matrix2image(active));

    figure(5);
    imshow(matrix2image(lowpass_img));
    title('The lowpass-filtered image');
    
    figure(6);
    title('The visual amplitude in the highpass-filtered frequency domain');
    hold on;
    imshow(uint8(highpass_amplitude));
    
    figure(7);
    title('The highpass-filter in the original Fourier frequency space');
    hold on;
    imshow(matrix2image(-active));

    figure(8);
    imshow(matrix2image(highpass_img));
    title('The highpass-filtered image');
end