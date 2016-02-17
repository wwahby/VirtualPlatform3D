function x = Noise_solver(C, Y, current, current_prev, x_prev, r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here    
% go to matrix build to find why the matrix is like these
    spparms('spumoni', 0);
    A = C/r+Y;
    B = C/r-Y;
    
    const = -(current + current_prev )/2;
    
    right = const + B*x_prev;
        
    tic;
    x = A\right;
    toc;
        
    %isequal(x, x_new)
end