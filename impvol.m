clc;
clear;
hold on;

steps = 300;
x = zeros(1, steps); 
y_cev = zeros(1, steps);
y_bls = zeros(1, steps);

sig = .25;                  % Stock volatility

for i = 1 : steps
    % Calculate initial call option price with strike K
    X   = 200;
    n   = 1400;
    m   = 140;
    K   = (i * X / steps);          % Strike price;
    T   = 1.0;                      % Years until expiry
    g   = @(x) max(0, x - K);       % Option function
    r   = 0.01;                     % Interest rate
    del = 1;                        % CEV Delta

    [u, time, space] = cn(T, X, n, m, g, r, del, sig);
    s = space(round(m / 10));
    
    cev_price = exp(-r * T) * u(n + 1, round(m / 10));
    bls_price = blsprice(s, K, r, T, sig);

    x(i) = K;
    y_cev(i) = blsimpv(s, K, r, T, cev_price);
    y_bls(i) = blsimpv(s, K, r, T, bls_price);
end

plot(x, y_cev, 'b');
plot(x, y_bls, 'k');

axis([5 (2 * s) 0 0.5]);