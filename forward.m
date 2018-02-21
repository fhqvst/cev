function [u, time, space] = forward(T, X, n, m, g, r, del, sig)
%FORWARD Summary of this function goes here
%   Detailed explanation goes here

dt = T/n;
dx = X/m;
d  = dt/dx^2;

if d < 0.5
    disp("Warning: d =" + num2str(d));
end

u = zeros(n + 1, m + 1);
time   = zeros(1, n + 1);
space  = zeros(1, m + 1);

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
    u(i, m + 1) = max(0, g(space(m + 1)));
    u(i, 1)     = max(0, g(space(1)));
end

for i = 1 : n
    for j = 2 : m
        x = space(j);
        a = r * x * (u(i, j + 1) - u(i, j - 1)) / (2 * dx);
        b = (sig^2 / 2) * x^(2 * del) * ((u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / dx^2);
        
        u(i + 1, j) = dt * (a + b) + u(i, j);
    end
end

end
