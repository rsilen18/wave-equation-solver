% Visualization of numerical Solution for BVP for 1-D Wave Equation:
%     Interval: x in (0,1)
%     u(x=0,t) = 0
%     u_x(x=1,t) = 0
%     Initial Value: u_0 = sin(3pi/2*x)
%     Initial Velocity: u_x0 = -sin(3pi/2*x)
% Author: Alexandre Ait-Ettajer
% Date: May 19, 2020 

N = 10;
M = 100;
t_f = 4;
a = 0;
b = 1;
dx = (b-a)/N;
x = linspace(a,b+dx,(N+2)); % Add Ghost Point
dt = t_f/M;
assert(x(2)-x(1) == dx);

xt = ones(N,M);

u_old = solver(;

for n = 1:M
    
end