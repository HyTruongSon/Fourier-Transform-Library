% Software: Discrete Fourier Transform 1D
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_dft_1d_operator()
    % Dimension
    N = 256;
    
    % Compute the Fourier operator
    [re1, im1] = operator_dft_1d(N); % Use this for further computation
    [re2, im2] = visual_operator_dft_1d(N); % Only use this for visualization, not computation
    
    % Visualization
    figure(1);
    subplot(1, 2, 1);
    imshow(matrix2image(re2));
    title('Real part - Cosine');
    subplot(1, 2, 2);
    imshow(matrix2image(im2));
    title('Imaginary part - Sine');
    
    % Compute the Fourier inverse operator
    [re3, im3] = operator_idft_1d(N); % Use this for further computation
    [re4, im4] = visual_operator_dft_1d(N); % Only use this for visualization, not computation
    
    % Visualization
    figure(2);
    subplot(1, 2, 1);
    imshow(matrix2image(re4));
    title('Real part - Cosine');
    subplot(1, 2, 2);
    imshow(matrix2image(im4));
    title('Imaginary part - Sine');
end