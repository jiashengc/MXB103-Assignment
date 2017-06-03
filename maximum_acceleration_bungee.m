function [a, maxacc] = maximum_acceleration_bungee(v,h,n)
%MAXIMUM ACCELERATION BUNGEE First order backward difference approximation for 
%acceleration for the bungee. [a, maxacc] = maximum_acceleration_bungee(v,h,n) returns 
%the maximum acceleration of the bungee jumper at all positions

a = zeros(length(n));

for j = 1:n
    a(j + 1) = (v(j + 1) - v(j))./h;
end

maxacc = max(abs(a));

if maxacc > 19.62
   warning('Maximum acceleration exceeds 2g/s')
end
