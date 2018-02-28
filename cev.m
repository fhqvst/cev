% Plot CEV prices and Black-Scholes prices for vanilla options.
%figure;
clc;
clear;
hold on;

plotStyleCall = {'--b','-b','-.b', 'vb', '~b'};
plotStylePut  = {'--r','-r','-.r', 'vr', '~r'};

T   = 1/2;
X   = 100;
n   = 100;
m   = 500;
K   = 10;
g   = @(x) x - K;
r   = 0.02;
deltas = {.75, 1, 1.1};
sig = 0.5;

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

% Plot Black-Scholes
[bsc, bsp] = blsprice(space, K, r, T, sig);

plot(x, bsc, '.k');
plot(x, bsp, '.k');

% Plot strike price
line([K K], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');
% line([(3/5 * K) (3/5 * K)], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');
% line([(7/5 * K) (7/5 * K)], [0 (K + 0.1)], 'Color', 'k', 'LineStyle', ':');

% Set axis
axis([0 (2 * K) 0 K]);
% axis auto;

hold off;

title('Vanilla put and call option pricing with the CEV model for strike K = 10');
legend('delta = 0.75', 'delta = 1.00', 'delta = 1.10', 'delta = 0.75', 'delta = 1.00', 'delta = 1.10', 'Black-Scholes');
xlabel('Stock price');
ylabel('Option price');


paper_width = 6;
paper_height = 6;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, paper_width, paper_height], 'PaperUnits', 'Inches', 'PaperSize', [paper_width, paper_height])
saveas(gcf,'figures/cev.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUT-CALL PARITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

del = 0.5;
[c_u, c_time, c_space] = cn(T, X, n, m, g, r, del, sig);
[p_u, p_time, p_space] = cn(T, X, n, m, @(x) -g(x), r, del, sig);

c_x = zeros(1, m + 1); p_x = zeros(1, m + 1);
c_y = zeros(1, m + 1); p_y = zeros(1, m + 1);

parity_errors = zeros(1, m + 1);

for j = 1 : m + 1
    c_x(j) = c_space(j);
    p_x(j) = p_space(j);

    c_y(j) = exp(-(r * T)) * c_u(n + 1, j);
    p_y(j) = exp(-(r * T)) * p_u(n + 1, j);
    
    error = abs(c_y(j) - p_y(j) - c_space(j) + K * exp(-r * T));
    parity_errors(j) = error;
end

avg_parity_error = mean(parity_errors);

plot(c_x, parity_errors);
hold on;
line([0 120], [avg_parity_error avg_parity_error], 'Color', 'k', 'LineStyle', ':');
hold off;

title('Error in the put-call parity for numerically calculated CEV model');
legend('error', 'average error', 'Location', 'northwest');
xlabel('Stock price');
ylabel('Error');
axis([0 110 0 1.1]);

paper_width = 6;
paper_height = 6;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, paper_width, paper_height], 'PaperUnits', 'Inches', 'PaperSize', [paper_width, paper_height])
saveas(gcf,'figures/parity.pdf')

disp(' ');
disp('Done!')