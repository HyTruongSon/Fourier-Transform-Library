% Software: Fast Fourier Transform 2D for Image Processing
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_2d_gaussian_denoising_by_lowpass_filter()
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    
    % Constants
    noise_magnitude = 20;
    radius = 100; % Counted by pixels. We only keep the frequency inside the radius from the (0, 0)
    
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
    
    % Compute the circle for Lowpass filter
    N = size(noised_img, 1);
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
    lowpass_amplitude = amp_n;    
    lowpass_amplitude(active == 0) = 0.0;
    
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
    lowpass_re = re_n;
    lowpass_im = im_n;
    
    lowpass_re(active == 0) = 0.0;
    lowpass_im(active == 0) = 0.0;
    
    % Reconstruct the image
    [lowpass_img] = ifft_2d(lowpass_re, lowpass_im);

    % Show the denoising results
    figure(5);
    title('The visual amplitude in the lowpass-filtered frequency domain');
    hold on;
    imshow(uint8(lowpass_amplitude));
    
    figure(6);
    title('The lowpass-filter in the Fourier frequency space');
    hold on;
    imshow(matrix2image(active));
    
    figure(7);
    imshow(matrix2image(lowpass_img));
    title('The lowpass-filtered image');
end