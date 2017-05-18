function [t, y, v, h] = rk4_bungee(T, n, g, C, K, L)


% Calculate subinterval width h
h = T / n;

% Create time array t
t = 0:h:T;

% Initialise solution arrays y and v
y = zeros(1,n+1);
v = zeros(1,n+1);

f = @(y,v) v;

% Newton's second law of motion
g = @(y,v) g - (C * abs(v) * v) - max(0, K * (y - L));

% Iterations
for i = 1:n
    k_1 = h * f(y(i), v(i));
    j_1 = h * g(y(i), v(i));
    
    k_2 = h * f(y(i) + k_1/2, v(i) + j_1/2);
    j_2 = h * g(y(i) + k_1/2, v(i) + j_1/2);
    
    k_3 = h * f(y(i) + k_2/2, v(i) + j_2/2);
    j_3 = h * g(y(i) + k_2/2, v(i) + j_2/2);
    
    k_4 = h * f(y(i) + k_3, v(i) + j_3);
    j_4 = h * g(y(i) + k_3, v(i) + j_3);
    
    y(i+1) = y(i) + (1/6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
    v(i+1) = v(i) + (1/6) * (j_1 + 2 * j_1 + 2 * j_3 + j_4);
end