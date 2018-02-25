function [u, time, space] = cn(T, X, n, m, g, r, del, sig)
%BACKWARD Summary of this function goes here
%   Detailed explanation goes here

dt = T/n;
dx = X/m;

d1  = dt/dx;
d2  = dt/(dx^2);

u = zeros(n + 1, m + 1);
time   = zeros(1, n + 1);
space  = zeros(1, m + 1);
A      = zeros(m + 1, m + 1);
B      = zeros(m + 1, m + 1);

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

% Values for A
for k = 1 : m + 1
    xj = space(k);
    
    a = r * xj * d1 / 4 - (sig^(2) * xj^(2 * del) * d2 / 4);
    b = 1 + (sig^(2) * xj^(2 * del) * d2 / 2);
    c = -r * xj * d1 / 4 - (sig^(2) * xj^(2 * del) * d2 / 4);

    if k > 1
        A(k, k - 1) = a;
    end

    A(k, k) = b;
    
    if k < m + 1
        A(k, k + 1) = c;
    end
end

% Values for B
for k = 1 : m + 1
    xj = space(k);

    a = (sig^(2) * xj^(2 * del) * d2 / 4) - (r * xj * d1 / 4);
    b = 1 - (sig^(2) * xj^(2 * del) * d2 / 2);
    c = (r * xj * d1 / 4) + (sig^(2) * xj^(2 * del) * d2 / 4);
    
    if k > 1
        B(k, k - 1) = a;
    end

    B(k, k) = b;
    
    if k < m + 1
        B(k, k + 1) = c;
    end
end

A_ = transpose(inv(A));
for i = 1 : n
    u(i + 1, :) = transpose(B * transpose(u(i, :))) * A_;
end

end
