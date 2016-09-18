% Software: Fast Fourier Transform 2D
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_2d_gaussian_denoising_by_smoothing()
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    
    % Constants
    noise_magnitude = 20;
    smoothing_radius = 5;
    smoothing_variance = 1.0;
    
    % The original image
    original_img = imread(img_fn);
    if (size(original_img, 3) == 3)
        original_img = rgb2gray(imread(img_fn));
    end
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    % FFT 2D of the original image
    [re_o, im_o] = fft_2d(double(original_img));
    [amp_o] = visual_amplitude_2d(re_o, im_o); % Only for visualization
    
    figure(2);
    title('The visual amplitude in the frequency domain (the original image)');
    hold on;
    imshow(uint8(amp_o));
    
    % The noised image
    noise = noise_magnitude * randn(size(original_img));
    noised_img = double(original_img) + noise;
    noised_img(noised_img < 0) = 0;
    noised_img(noised_img > 255) = 255;
    
    figure(3);
    imshow(uint8(noised_img));
    title(['The noised image with the gaussian noise (magnitude = ', int2str(noise_magnitude), ')']);
    
    % FFT 2D of the noised image
    [re_n, im_n] = fft_2d(noised_img);
    [amp_n] = visual_amplitude_2d(re_n, im_n); % Only for visualization
    
    figure(4);
    title('The visual amplitude in the frequency domain (the noised image)');
    hold on;
    imshow(uint8(amp_n));
    
    % The Gaussian matrix with the smoothing radius and variance
    gm = gaussian_matrix([smoothing_radius, smoothing_variance]);
    denoised_img = image_filter(double(noised_img), gm);
    
    figure(5);
    imshow(uint8(denoised_img));
    title(['The denoised image by Guassian matrix (radius = ', int2str(smoothing_radius), ', variance = ', int2str(smoothing_variance), ')']);

    % FFT 2D of the denoised image
    [re_d, im_d] = fft_2d(denoised_img);
    [amp_d] = visual_amplitude_2d(re_d, im_d); % Only for visualization
    
    figure(6);
    title('The visual amplitude in the frequency domain (the denoised image)');
    hold on;
    imshow(uint8(amp_d));
end