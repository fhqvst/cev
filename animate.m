% Animation how the numerical method for CEV approximates Black-Scholes.
clc;
figure

T   = 1;
X   = 2;
n   = 1000;
m   = 100;
K   = 0.5;
g   = @(x) x - K;
r   = 0.01;
del = 1;
sig = 0.25;

disp('Running...');
disp('');

v = 0.1;
step = round(1 + n * v / 4);

for i = 1 : step : n + 1

    [u, time, space] = cn(T, X, n, m, g, r, del, sig);

    S = 1 : m + 1;
    x = zeros(m + 1);
    y = zeros(m + 1);

    for j = 1 : m + 1
        x(j) = space(j);
        y(j) = exp(-(r * time(i))) * u(n + 2 - i, j);
    end

    plot(x, y, 'k');
    
    for j = 1 : m + 1
        [call, put] = blsprice(space(j), K, r, time(n + 2 - i), sig);
        y(j) = call;
    end
    hold on;
    plot(x, y, ':k');
    hold off;
    
    axis auto;
    drawnow;

end

disp('Done!')