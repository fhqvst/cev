clc;
clear;
% figure;
hold on;

plotStyleCall = {'-.b','.b','--b', 'vb', '~k'};
plotStylePut = {'-.r','.r','--r', 'vr', '~k'};

T   = 1;
X   = 1;
n   = 10000;
m   = 100;
K   = 0.5;
g   = @(x) x - K;
r   = 0.01;
deltas = {0.5, 1, 1.5};
sig = 0.25;

disp('Running...');
disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL (CEV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot CEV
disp('Calls:');
disp('');
for k = 1 : length(deltas)
    del = deltas{k};

    [u, time, space] = cn(T, X, n, m, g, r, del, sig);

    S = 1 : m + 1;
    x = zeros(m + 1);
    y = zeros(m + 1);

    for j = 1 : m + 1
        x(j) = space(j);
        y(j) = exp(-(r * time(1))) * u(n + 1, j);
    end
    
    plot(x, y, plotStyleCall{k});
    disp(num2str(del) + " = " + num2str(plotStyleCall{k}));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT (CEV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Puts:');
disp('');

g   = @(x) K - x;
for k = 1 : length(deltas)
    del = deltas{k};

    [u, time, space] = cn(T, X, n, m, g, r, del, sig);

    S = 1 : m + 1;
    x = zeros(m + 1);
    y = zeros(m + 1);

    for j = 1 : m + 1
        x(j) = space(j);
        y(j) = exp(-(r * time(1))) * u(n + 1, j);
    end
    
    plot(x, y, plotStylePut{k});
    disp(num2str(del) + " = " + num2str(plotStylePut{k}));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot strike price
disp("K = " + num2str(K));
line([K K], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');

% Plot Black-Scholes
bsc = zeros(m + 1);
bsp = zeros(m + 1);
for j = 1 : m + 1
    [call, put] = blsprice(space(j), K, r, time(n + 1), sig);
    bsc(j) = call;
    bsp(j) = put; 
end

plot(x, bsc, '-b');
plot(x, bsp, '-r');

axis([0 1 0 0.6]);
hold off;

disp('Done!')