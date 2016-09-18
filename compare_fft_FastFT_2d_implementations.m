% Software: Fast Fourier Transform 2D
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = compare_fft_FastFT_2d_implementations()
	fprintf('Comparison the performance of fft and FastFT 2D:');
    
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    original_img = imread(img_fn);
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    if (size(original_img, 3) == 3)
        original_img = rgb2gray(imread(img_fn));
    end
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    fprintf('- fft_2d and ifft_2d\n');
    fprintf('  Time complexity: O(N^2logN)\n');
    fprintf('  Space complexity: O(N^2)\n');
    
    [re1, im1] = fft_2d(double(original_img));
    [img1] = ifft_2d(re1, im1);
    
    figure(2);
    imshow(uint8(img1));
    title('The reconstructed image from fft-2d.cpp and ifft-2d.cpp');
    
    fprintf('- FastFT_2d and iFastFT_2d\n');
    fprintf('  Time complexity: O(N^2logN)\n');
    fprintf('  Space complexity: O(N^2logN)\n');
    
    [re2, im2] = FastFT_2d(double(original_img));
    [img2] = iFastFT_2d(re2, im2);
    
    figure(3);
    imshow(uint8(img2));
    title('The reconstructed image from FastFT-2d.cpp and iFastFT-2d.cpp');
end