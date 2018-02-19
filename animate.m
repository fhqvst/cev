clc;
figure

T   = 1;
X   = 1;
n   = 1000;
m   = 30;
K   = 0.5;
g   = @(x) x - K;
r   = 0.01;
del = 0.5;
sig = 0.25;

disp('Running...');
disp('');

v = 0.1;
step = round(1 + n * v / 10);

for i = 1 : step : n + 1

    [u, time, space] = forward(T, X, n, m, g, r, del, sig);

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
    
    axis([0 1 0 .6]);
    drawnow;

end

disp('Done!')