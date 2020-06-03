% driver code for waves_fdm_1d.m
N = 100;    % resolution in x-dimension
c = 1;      % wave speed
t_f = 10;   % final time
sigma = 0.5; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
order = 4;  % 2nd or 4th order precision
plot_flag = true;    % flag for turning on plot animation

% Original problem
[u,e] = waves_fdm_1d(N,sigma,c,t_f,icase,plot_flag,order);

hold on
x = linspace(0,3*pi/2,N+1);
plot(x,waves_analytic_1d(@(x) sin(x),@(x) -sin(x),c,x,t_f));
hold off

% Convergence study
icase = 1;
plot_flag = false;
h = [.1 .01 .001 ];
N = [10 100 1000 ];
errors = zeros(1,3);
for i = 1:3
    disp(i);
    [u,e] = waves_fdm_1d(N(i),sigma,c,t_f,icase,plot_flag,order);
    errors(i) = max(e);
end
h2 = h.^2;
h4 = h.^4;
figure(2)
loglog(h,errors,'o',h,h2,h,h4);
xlabel("h");
ylabel("|e|_{\infty}");
title("error vs. dx");