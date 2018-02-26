function [u, time, space] = backward(T, X, n, m, g, r, del, sig)
%BACKWARD BACKWARD in time, centered in space.
%   Backward euler method in the time derivative, centered in the space
%   derivative.

dt = T/n;
dx = X/m;

d1  = dt/dx;
d2  = dt/(dx^2);

u = zeros(n + 1, m + 1);
time   = zeros(1, n + 1);
space  = zeros(1, m + 1);
A      = zeros(m + 1, m + 1);

for i = 2 : n + 1
    time(i) = time(i - 1) + dt;
end

for j = 2 : m + 1
    space(j) = space(j - 1) + dx;
end

% Value for initial time
for j = 1 : m + 1
    u(1, j) = max(0, g(space(j)));
end

% Value for initial space and terminal space
for i = 1 : n + 1
    u(i, 1)     = max(0, g(space(1) / exp(-r * time(i)) * exp(-r * time(i))));
    u(i, m + 1) = max(0, g(space(m + 1) / exp(-r * time(i)) * exp(-r * time(i))));
end

% Fix A values for given boundaries
A(1, 1) = 1; 
A(m + 1, m + 1) = 1;

% Values for A
for k = 2 : m
    xj = space(k);

    a = r * xj * d1 / 2 - (sig^(2) * xj^(2 * del) * d2 / 2);
    b = 1 + (sig^(2) * xj^(2 * del) * d2);
    c = -(r * xj * d1 / 2 + sig^(2) * xj^(2 * del) * d2 / 2);

    A(k, k - 1) = a;
    A(k, k) = b;
    A(k, k + 1) = c;
end

A_ = transpose(inv(A));
for i = 1 : n
    u(i + 1, :) = u(i, :) * A_;
end

end
