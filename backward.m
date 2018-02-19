function [u, time, space] = backward(T, X, n, m, g, r, del, sig)
%BACKWARD Summary of this function goes here
%   Detailed explanation goes here

[u, time, space] = forward(T, X, n, m, g, r, del, sig);

end