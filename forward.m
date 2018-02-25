function [u, time, space] = forward(T, X, n, m, g, r, del, sig)
%FORWARD Summary of this function goes here
%   Detailed explanation goes here

dt = T/n;
dx = X/m;

d1 = dt/dx;
d2 = dt/(dx^2);

disp(['d1: ', num2str(d1)]);
disp(['d2: ', num2str(d2)]);

if (d1 > 0.00068)
    disp(['Warning: d1 > 0.0068 (', num2str(d1), ')']);
end
if (d2 > 0.0004)
    disp(['Warning: d2 > 0.00046923 (', num2str(d2), ')']);
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
    u(i, 1)     = max(0, g(space(1) / exp(-r * time(i)) * exp(-r * time(i))));
    u(i, m + 1) = max(0, g(space(m + 1) / exp(-r * time(i)) * exp(-r * time(i))));
end

for i = 1 : n
    for j = 2 : m
        xj = space(j);
        
        a = r * xj * d1 * (u(i, j + 1) - u(i, j - 1)) / 2;
        b = sig^(2) * xj^(2 * del) * d2 * (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / 2;
        
        u(i + 1, j) = u(i, j) + a + b;
    end
end

end
