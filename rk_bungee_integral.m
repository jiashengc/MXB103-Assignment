function [t, y, v, h, x] = rk_bungee_integral(T, n, g, C, K, L)
% rk_bungee Runge-Kutta method for the bungee jumping model
% [t, y, v, h] = rk_bungee(T, n, g, C, K, L) performs Euler's method on
% the bungee jumping model, taking n steps from t = 0 to t = T.
% The initial conditions are y(0) = 0 and v(0) = 0.
% The inputs g, C, K and L are parameters from the model 
% (see project description). The outputs are the time array t, the solution
% arrays y and v, and the subinterval width h. To be complimented by a 
% time vs displacement plot.

    % Calculate subinterval width h
    h = T / n;

    % Create time array t
    t = 0:h:T;
        
    % Initialise solution arrays y and v
    y = zeros(1,n+1);
    v = zeros(1,n+1);
    x = zeros(1,n+1);

    f = @(t,y,v) v;
    
    g = @(t,y,v) g - C*abs(v)*v - max(0, K*(y - L));
    
    % Perform iterations
    for j = 1:n
        k1 = h * f(t(j), y(j), v(j));
        n1 = h * g(t(j), y(j), v(j));
        
        k2 = h * f(t(j) + h/2, y(j) + k1/2, v(j) + n1/2);
        n2 = h * g(t(j) + h/2, y(j) + k1/2, v(j) + n1/2);
        
        k3 = h * f(t(j) + h/2, y(j) + k2/2, v(j) + n2/2);
        n3 = h * g(t(j) + h/2, y(j) + k2/2, v(j) + n2/2);
        
        k4 = h * f(t(j) + h, y(j) + k3, v(j) + n3);
        n4 = h * g(t(j) + h, y(j) + k3, v(j) + n3);
        
       y(j+1) = y(j) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
       v(j+1) = v(j) + (1/6)*(n1 + 2*n2 + 2*n3 + n4);
       x(j+1) = x(j) + abs((1/6)*(k1 + 2*k2 + 2*k3 +k4));
    end
          
   