% Software: Fast Fourier Transform 2D for Image Processing
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_fft_2d_hybrid_frequecy_image()
    % Image file names
    img_fn_1 = 'lena_gray.bmp'; % Size 512 x 512
    img_fn_2 = 'girlface.bmp'; % Size 512 x 512
    
    % Constants
    radius = 50;
    
    % Convert the images into gray-scale
    img_1 = imread(img_fn_1);
    if size(img_1, 3) > 1
        img_1 = rgb2gray(img_1);
    end
    
    img_2 = imread(img_fn_2);
    if size(img_2, 3) > 1
        img_2 = rgb2gray(img_2);
    end
    
    if size(img_1) ~= size(img_2)
        fprintf('The sizes of two images are not the same!');
        return;
    end
    
    % Show the first image
    figure(1);
    imshow(img_1);
    title('The first image');
    
    % FFT 2D of the first image
    [re_1, im_1] = fft_2d(double(img_1));
    [amp_1] = visual_amplitude_2d(re_1, im_1); % Only for visualization
    
    % Show the first image's frequency space
    figure(2);
    imshow(uint8(amp_1));
    title('The frequency space of the first image');
    
    % Show the first image
    figure(3);
    imshow(img_2);
    title('The second image');
    
    % FFT 2D of the second image
    [re_2, im_2] = fft_2d(double(img_2));
    [amp_2] = visual_amplitude_2d(re_2, im_2); % Only for visualization
    
    % Show the second image's frequency space
    figure(4);
    imshow(uint8(amp_2));
    title('The frequency space of the second image');
    
    % Compute the circle for Lowpass filter
    N = size(img_1, 1);
    active_amp = zeros(N, N);
    x0 = N / 2;
    y0 = N / 2;
    for x = x0 - radius : x0 + radius
        delta = radius * radius - (x - x0) * (x - x0);
        if delta >= 0
            r = int32(sqrt(delta));
            active_amp(x, y0 - r : y0 + r) = 1;
        end
    end
    
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
    
    % Filtering: The frequency of the first image is inside the circle,
    % outside is the second image's
    amp_12 = amp_1;
    amp_12(active_amp == 0) = amp_2(active_amp == 0);
    re_12 = re_1;
    re_12(active == 0) = re_2(active == 0);
    im_12 = im_1;
    im_12(active == 0) = im_2(active == 0);
    
    % Show the frequency of the 1-2 case
    figure(5);
    imshow(uint8(amp_12));
    title('The first image frequency is inside, outside is the second');
    
    % Inverse Fourier Transform
    [img_12] = ifft_2d(re_12, im_12);
    
    % Show the image 12
    figure(6);
    imshow(uint8(img_12));
    title('The first image frequency is inside, outside is the second');
    
    % Filtering: The frequency of the second image is inside the circle,
    % outside is the first image's
    amp_21 = amp_2;
    amp_21(active_amp == 0) = amp_1(active_amp == 0);
    re_21 = re_2;
    re_21(active == 0) = re_1(active == 0);
    im_21 = im_2;
    im_21(active == 0) = im_1(active == 0);
    
    % Show the frequency of the 1-2 case
    figure(7);
    imshow(uint8(amp_21));
    title('The second image frequency is inside, outside is the first');
    
    % Inverse Fourier Transform
    [img_21] = ifft_2d(re_21, im_21);
    
    % Show the image 21
    figure(8);
    imshow(uint8(img_21));
    title('The second image frequency is inside, outside is the first');
end