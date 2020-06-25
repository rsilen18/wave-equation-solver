% driver code for waves_fdm_1d.m
sigma = 0.5; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
order = 6;  % 2nd or 4th order precision
plot_flag = true;    % flag for turning on plot animation

% problem definition
[def.a,def.b,def.c,def.N,def.t_f,def.f,def.g,def.left,def.right] = waves_fdm_1d_defs(icase);

% Original problem
[u,e] = waves_fdm_1d(def,sigma,plot_flag,order);

% hold on
% x = linspace(0,3*pi/2,def.N+1);
% plot(x,waves_analytic_1d(@(x) sin(x),@(x) -sin(x),c,x,def.t_f));
% hold off
% plot(x,waves_tz_1d(def,sigma,true));

% Convergence study
icase = 1;
plot_flag = false;
h = [.1 .01 .001 ];
N = [10 100 1000 ];
errors = zeros(1,3);
for i = 1:3
    disp(i);
    def.N = N(i);
    [u,e] = waves_fdm_1d(def,sigma,plot_flag,order);
    errors(i) = max(e);
end
h2 = h.^2;
h3 = h.^3;
h4 = h.^4;
h5 = h.^5;
h6 = h.^6;
figure(2)
loglog(h,errors,'o',h,h2,h,h3,h,h4,h,h5,h,h6);
xlabel("h");
ylabel("|e|_{\infty}");
title("Truncation Error");
legend("error","h^2","h^3","h^4","h^5","h^6");