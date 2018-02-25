clc;
clear;
hold on;

K_max = 40;
x = zeros(1, K_max); 

y_cev_c = zeros(1, K_max);
y_cev_p = zeros(1, K_max);
y_bls = zeros(1, K_max);

for i = 1 : K_max
    % Calculate initial call option price with strike K
    K   = i;
    X   = 200;
    n   = 1000;
    m   = 600;
    T   = 1;                        % Years until expiry
    g   = @(x) max(0, x - K);       % Option function
    r   = 0.01;                     % Interest rate
    del = 1;                        % CEV Delta
    sig = .25;                      % Volatility

    [u_c, time_c, space_c] = cn(T, X, n, m, g, r, del, sig);

    si = round(m / 10);
    s0 = space_c(si);
    
    cev_price_c = exp(-r * T) * u_c(n + 1, si);
    bls_price = blsprice(s0, K, r, T, sig);

    x(i) = K;

    y_cev_c(i) = blsimpv(s0, K, r, T, cev_price_c, 1, 0, [], true);
    y_bls(i) = blsimpv(s0, K, r, T, bls_price);
end

plot(x, y_cev_c, 'b');
plot(x, y_bls, 'k');

axis auto;