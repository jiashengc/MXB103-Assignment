function [t, y, v, h, x] = rk_bungee_integral(T, n, g, C, K, L)
% rk_bungee Runge-Kutta method for the bungee jumping model
% [t, y, v, h] = rk_bungee(T, n, g, C, K, L) implements Euler's method on
% the bungee jumping model, taking n steps from t = 0 to t = T.
% The initial conditions are y(0) = 0 and v(0) = 0.
% The inputs g, C, K and L are parameters from the model 
% The outputs are the time array t, the solution
% arrays y and v, and the subinterval width h. To be complimented by a 
% time vs displacement plot.

    % Subinterval width h
    h = T / n;

    % Time array t
    t = 0:h:T;
        
    % Solution arrays y and v
    y = zeros(1,n+1);
    v = zeros(1,n+1);
    x = zeros(1,n+1);

    f = @(t,y,v) v;
    
    g = @(t,y,v) g - C*abs(v)*v - max(0, K*(y - L));
    
    % Iterations
    for i = 1:n
        k_1 = h * f(t(i), y(i), v(i));
        n_1 = h * g(t(i), y(i), v(i));
        
        k_2 = h * f(t(i) + h/2, y(i) + k_1/2, v(i) + n_1/2);
        n_2 = h * g(t(i) + h/2, y(i) + k_1/2, v(i) + n_1/2);
        
        k_3 = h * f(t(i) + h/2, y(i) + k_2/2, v(i) + n_2/2);
        n_3 = h * g(t(i) + h/2, y(i) + k_2/2, v(i) + n_2/2);
        
        k_4 = h * f(t(i) + h, y(i) + k_3, v(i) + n_3);
        n_4 = h * g(t(i) + h, y(i) + k_3, v(i) + n_3);
        
       y(i+1) = y(i) + (1/6)*(k_1 + 2*k_2 + 2*k3 + k_4);
       v(i+1) = v(i) + (1/6)*(n_1 + 2*n_2 + 2*n3 + n_4);
       x(i+1) = x(i) + abs((1/6)*(k_1 + 2*k_2 + 2*k3 +k_4));
    end
          
   