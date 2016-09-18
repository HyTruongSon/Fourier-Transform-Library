% Software: Discrete Fourier Transform 2D (for real signals only)
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_dft_2d()
    img_fn = 'index.jpeg'; % Size: 204 x 204
    
    original_img = imread(img_fn);
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    if (size(original_img, 3) == 3)
        original_img = rgb2gray(imread(img_fn));
    end
    
    % Discrete Fourier Transform 2D
    [re, im] = dft_2d(double(original_img));
    
    % Inverse Discrete Fourier Transform 2D
    [reversed_img] = idft_2d(re, im);
    
    figure(2);
    title('The reversed image after DFT 2D and Inverse DFT 2D');
    hold on;
    reversed_img = uint8(reversed_img);
    imshow(reversed_img);
    
    % Compute amplitude
    [original_amplitude] = amplitude_2d(re, im); % Correct amplitude
    [visual_amplitude] = visual_amplitude_2d(re, im); % Only for visualization
    
    figure(3);
    title('The visual amplitude in the frequency domain');
    hold on;
    imshow(uint8(visual_amplitude));
    
    % Compute phase
    [original_phase] = phase_2d(re, im); % Correct phase
    [visual_phase] = visual_phase_2d(re, im); % Only for visualization
    
    figure(4);
    title('The visual phase in the frequency domain');
    hold on;
    imshow(uint8(visual_phase));
end