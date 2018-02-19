function [u, time, space] = cn(T, X, n, m, g, r, del, sig)
%CN Summary of this function goes here
%   Detailed explanation goes here

% Run forward
[u_f, time_f, space_f] = forward(T, X, n, m, g, r, del, sig);

% Run backward
[u_b, time_b, space_b] = backward(T, X, n, m, g, r, del, sig);

% Crank-Nicholson
u = 1/2 * (u_f + u_b);
time = 1/2 * (time_f + time_b);
space = 1/2 * (space_f + space_b);

end

