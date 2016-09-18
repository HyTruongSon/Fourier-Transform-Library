% Software: Fast Fourier Transform 2D for Image Processing
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [img] = matrix2image(img)
    MAX = max(max(img));
    MIN = min(min(img));
    img = (img - MIN) / (MAX - MIN) * 255.0;
    img = uint8(img);
end