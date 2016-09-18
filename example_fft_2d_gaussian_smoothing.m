% Software: Fast Fourier Transform 2D for Image Processing
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_2d_gaussian_smoothing()
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    
    % Constants
    radius = 5; % Pixels
    variance = 1.0;
    
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
    
    % Gaussian smoothing
    filtered_img = image_filter(double(original_img), gaussian_matrix([radius, variance]));
    
    figure(3);
    imshow(uint8(filtered_img));
    title('The filtered image');
    
    % Fast Fourier Transform 2D of the filtered image
    [re, im] = fft_2d(double(filtered_img));
    
    % Compute amplitude of the filtered image
    [visual_amplitude] = visual_amplitude_2d(re, im); % Only for visualization
    
    figure(4);
    title('The visual amplitude in the frequency domain (filtered image)');
    hold on;
    imshow(uint8(visual_amplitude));
end