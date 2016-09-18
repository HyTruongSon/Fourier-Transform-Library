% Software: Discrete Cosine Transform 2D for JPEG
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_dct_2d_jpeg()
    % img_fn = 'lena.png'; % Size: 256 x 256
    img_fn = 'lena_gray.bmp'; % Size: 512 x 512
    
    original_img = imread(img_fn);
    
    figure(1);
    imshow(original_img);
    title('The original image');
    
    if (size(original_img, 3) == 3)
        original_img = rgb2gray(imread(img_fn));
    end
    
    % Discrete Cosine Transform 2D
    nRows = size(original_img, 1);
    nCols = size(original_img, 2);
    
    [dct] = dct_2d_jpeg(double(original_img), nRows, nCols);
    
    figure(2);
    imshow(matrix2image(dct));
    title('The Discrete Cosine Transform 2D - block size 8x8');
    
    % Quantization - JPEG Standard
    [luminance] = quantization_table(1);
    [chrominance] = quantization_table(2);
    
    [lu_dct] = quantize(dct, luminance);
    [ch_dct] = quantize(dct, chrominance);
    
    % +---------------------------------------------+
    % | Entropy coding - Human code - Lossless code |
    % +---------------------------------------------+
    
    % Inverse Quantization
    [lu_dct] = iquantize(lu_dct, luminance);
    [ch_dct] = iquantize(ch_dct, chrominance);
    
    % Inverse Discrete Cosine Transform 2D
    [full_img] = idct_2d_jpeg(dct, nRows, nCols);
    [lu_img] = idct_2d_jpeg(lu_dct, nRows, nCols);
    [ch_img] = idct_2d_jpeg(ch_dct, nRows, nCols);
    
    figure(3);
    % subplot(1, 3, 1);
    imshow(uint8(full_img));
    title('The reversed full image from the Inverse Discrete Cosine Transform 2D');
    
    figure(4);
    % subplot(1, 3, 2);
    imshow(uint8(lu_img));
    title('The JPEG image with Quantization Luminance');
    
    figure(5);
    % subplot(1, 3, 3);
    imshow(uint8(ch_img));
    title('The JPEG image with Quantization Chrominance');
    
    figure(6);
    imshow(uint8(abs(double(original_img) - lu_img)));
    title('Difference between the original image and the Luminance quantized image');
    
    figure(7);
    imshow(uint8(abs(double(original_img) - ch_img)));
    title('Difference between the original image and the Chrominance quantized image');
end