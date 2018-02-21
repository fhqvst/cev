clc;
clear;
% figure;
hold on;

plotStyleCall = {'-.b','.b','--b', 'vb', '~b'};
plotStylePut  = {'-.r','.r','--r', 'vr', '~r'};

T   = 1;
X   = 1;
n   = 15000;
m   = 400;
K   = 0.5;
g   = @(x) x - K;
r   = 0.01;
deltas = {0.5, 1.0, 1.5};
sig = 0.25;

disp('Running...');
disp(' ');
disp(['Strike: ', num2str(K)]);
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL (CEV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Calls: (deltas)');

for k = 1 : length(deltas)
    del = deltas{k};

    [u, time, space] = cn(T, X, n, m, g, r, del, sig);

    x = zeros(1, m + 1);
    y = zeros(1, m + 1);

    for j = 1 : m + 1
        x(j) = space(j);
        y(j) = exp(-(r * T)) * u(n + 1, j);
    end
    
    plot(x, y, plotStyleCall{k});
    disp("  " + num2str(del) + " = " + plotStyleCall{k});
end

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT (CEV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Puts: (deltas)');

for k = 1 : length(deltas)
    del = deltas{k};

    [u, time, space] = cn(T, X, n, m, @(x) -g(x), r, del, sig);

    x = zeros(1, m + 1);
    y = zeros(1, m + 1);

    for j = 1 : m + 1
        x(j) = space(j);
        y(j) = exp(-(r * T)) * u(n + 1, j);
    end
    
    plot(x, y, plotStylePut{k});
    disp("  " + num2str(del) + " = " + plotStylePut{k});
end

disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot strike price
line([K K], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');
line([(3/5 * K) (3/5 * K)], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');
line([(7/5 * K) (7/5 * K)], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');

% Plot Black-Scholes
bsc = zeros(1, m + 1);
bsp = zeros(1, m + 1);
for j = 1 : m + 1
    [call, put] = blsprice(space(j), K, r, T, sig);
    bsc(j) = call;
    bsp(j) = put; 
end

plot(x, bsc, '-k');
plot(x, bsp, '-k');

axis([0 1 0 (K + 0.1)]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT-CALL PARITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

del = 0.5;
[c_u, c_time, c_space] = cn(T, X, n, m, g, r, del, sig);
[p_u, p_time, p_space] = cn(T, X, n, m, @(x) -g(x), r, del, sig);

c_x = zeros(1, m + 1); p_x = zeros(1, m + 1); 
c_y = zeros(1, m + 1); p_y = zeros(1, m + 1);

parity = 0;

for j = 1 : m + 1
    c_x(j) = c_space(j);
    p_x(j) = p_space(j);

    c_y(j) = exp(-(r * T)) * c_u(n + 1, j);
    p_y(j) = exp(-(r * T)) * p_u(n + 1, j);
    
    parity = parity + abs(c_y(j) - p_y(j) - c_space(j) + K * exp(-r * T));
end

disp(['Parity error = ', num2str(parity / (m + 1))]);

disp(' ');
disp('Done!')