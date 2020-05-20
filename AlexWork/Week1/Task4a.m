% Visualization of numerical Solution for BVP for 1-D Wave Equation:
%     Interval: x in (0,1)
%     u(x=0,t) = 0
%     u_x(x=1,t) = 0
%     Initial Value: u_0 = sin(3pi/2*x)
%     Initial Velocity: u_x0 = -sin(3pi/2*x)
% Author: Alexandre Ait-Ettajer
% Date: May 19, 2020 

% Initialize Values for Particular solve
N = 100;
M = 400;
t_f = 4;
c = .5;
a = 0;
b = 1;
dx = (b-a)/N;
x = linspace(a,b+dx,(N+2)); % Add Ghost Point
t = linspace(0,t_f,M);
dt = t_f/(M-1);
assert(t(2)-t(1) == dt);
assert(x(2)-x(1) == dx);
assert(c^2*dt^2/dx^2 < 1);

xt = ones(M,N+2);

for n = 1:M
    u_n = 0;
    u_n_1 = 0;
    if (n == 2)
        u_n = xt(n-1,:);
        u_n_1 = xt(n-1,:);
    elseif (n ~= 1)
        u_n_1 = xt(n-2,:);
        u_n = xt(n-1,:);
    end
    xt(n,:) = solver(u_n,u_n_1,dt*(n-1),dt,x,dx,c);
    figure(1)
    plot(x,xt(n,:));
    xlim([a b])
    ylim([-2 2])
    title('Wave Equation Solution: BVP Snapshot')
    xlabel('x')
    ylabel('u')
    pause(0.01) 
end

mesh(x,t,xt)
xlim([a b])
title('Wave Equation Solution: BVP Over Time')
xlabel('x (Distance)')
ylabel('t (Time)')
zlabel('u (Function Value)')
