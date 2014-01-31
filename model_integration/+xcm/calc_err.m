function [err_raw err_norm] = calc_err(a,b)
% calculate the error between two vectors
% a and b must be 1D vectors, but do not need to have the same lengths

% Pad shorter vector with zeros
ldiff1 = length(a) - length(b);
ldiff2 = length(b) - length(a);
a = [a zeros(1,ldiff2)];
b = [b zeros(1,ldiff1)];

err_raw = abs(b-a);
err_norm = err_raw./a;