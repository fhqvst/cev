% Plot Black-Scholes implied volatility, where price is given by the CEV-model.
% clc;
clear;
figure;
hold on;

plotStyle = {'--b','-b','-.b', 'vb', '~b'};


% deltas = {1.499, 1.5, 1.501}; % Helt rak, startar på omkring 1.39
deltas = {0.999, 1, 1.001}; % Böjen uppåt, startar på omkring 0.494
% deltas = {0.499, 0.5, 0.501}; % Böjen nedåt, startar på omkring 0.053

ddeltas = {
    {1.499, 1.5, 1.501}, % Helt rak, startar på omkring 1.39
    {0.999, 1, 1.001},   % Böjen uppåt, startar på omkring 0.494
    {0.499, 0.5, 0.501} % Böjen nedåt, startar på omkring 0.053
};

for h = 1 : length(deltas)
    deltas = ddeltas{h};
    subplot(1, length(ddeltas), h);
    hold on;

    for k = 1 : length(deltas)
        K_max = 120;
        K_min = 80;

        x = ones(1, K_max - K_min); 

        y_cev_c = ones(1, K_max - K_min);
        y_cev_p = ones(1, K_max - K_min);
        y_bls = ones(1, K_max - K_min);

        for i = K_min : K_max
            K   = i;
            X   = 200;
            n   = 100;
            m   = 500;
            T   = 1;                        % Years until expiry
            g   = @(x) max(0, x - K);       % Option function
            r   = 0.01;                     % Interest rate
            del = deltas{k};                % CEV Delta
            sig = .25;                      % Volatility

            [u_c, time_c, space_c] = cn(T, X, n, m, g, r, del, sig);

            si = round(m / 2);
            s0 = space_c(si);

            cev_price_c = exp(-r * T) * u_c(n + 1, si);
            bls_price = blsprice(s0, K, r, T, sig);

            x(i - K_min + 1) = K;

            y_cev_c(i - K_min + 1) = blsimpv(s0, K, r, T, cev_price_c);
            y_bls(i - K_min + 1)   = blsimpv(s0, K, r, T, bls_price);
        end
        disp(['Price: ', num2str(s0)]);

        plot(x, y_cev_c);
        %line([s0 s0], [0 .1], 'Color', 'k', 'LineStyle', ':');

        xlabel('Strike price');
        ylabel('Implied volatility');
        
        if h == 2
            title('Black-Scholes implied volatility for vanilla call option with price given by the CEV model');
        end

        axis auto;
    end
    legend("delta = " + deltas{1}, "delta = " + deltas{2}, "delta = " + deltas{3});
end

paper_width = 24;
paper_height = 6;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, paper_width, paper_height], 'PaperUnits', 'Inches', 'PaperSize', [paper_width, paper_height])
saveas(gcf,'figures/implied-volatility.pdf')

hold off;

%plot(x, y_bls, 'k');