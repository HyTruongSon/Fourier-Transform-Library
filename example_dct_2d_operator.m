% Software: Discrete Cosine Transform Operator 2D for JPEG
% Author: Hy Truong Son
% Position: PhD Student
% Institution: Department of Computer Science, The University of Chicago
% Email: sonpascal93@gmail.com, hytruongson@uchicago.edu
% Website: http://people.inf.elte.hu/hytruongson/
% Copyright 2016 (c) Hy Truong Son. All rights reserved.

function [] = example_dct_2d_operator()
    % Number of rows
    M = 8;
    
    % Number of columns
    N = 8;
    
    % Discrete Cosine Transform Operator 2D
    [dct_op] = operator_dct_2d(M, N);
    
    figure(1);
    title('DCT operator 2D');
    for x = 1 : M
        for y = 1 : N
            p = (x - 1) * M + y;
            subplot(M, N, p);
            block = dct_op((x - 1) * M + 1 : x * M, (y - 1) * N + 1 : y * N);
            if p > 1
                imshow(matrix2image(block));
            else
                block(:, :) = 255;
                imshow(block);
            end
            drawnow;
        end
    end
    
    % Inverse Discrete Cosine Transform Operator 2D
    [idct_op] = operator_idct_2d(M, N);
    
    figure(2);
    title('IDCT operator 2D');
    for x = 1 : M
        for y = 1 : N
            p = (x - 1) * M + y;
            subplot(M, N, p);
            block = idct_op((x - 1) * M + 1 : x * M, (y - 1) * N + 1 : y * N);
            imshow(matrix2image(block));
            drawnow;
        end
    end
    
%     fprintf('Difference: %.6f\n', sum(sum(abs(dct_op - idct_op))));
end