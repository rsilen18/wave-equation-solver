% driver code for waves_fdm_1d.m
sigma = 0.2; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
order = 6;  % 2nd or 4th order precision
plot_flag = true;    % flag for turning on plot animation
tz_flag = false;     % flag for twilight zone calculations

% problem definition
[def.a,def.b,def.c,def.N,def.t_f,def.f,def.g,def.left,def.right] = waves_fdm_1d_defs(icase);

% Original problem
% [u,e] = waves_fdm_1d(def,sigma,plot_flag,order,tz_flag);

% hold on
% x = linspace(0,3*pi/2,def.N+1);
% plot(x,waves_analytic_1d(@(x) sin(x),@(x) -sin(x),c,x,def.t_f));
% hold off
% plot(x,waves_tz_1d(def,sigma,true));



% Convergence study
icase = 1;
plot_flag = false;
h = [.1 .02 .01 .002 .001 ];
N = 1./h;
errors = zeros(1,size(h,2));
errors_tz = zeros(1,size(h,2));
for i = 1:size(h,2)
    disp(i);
    def.N = N(i);
    [u,e] = waves_fdm_1d(def,sigma,plot_flag,order,false);
    [u_tz,e_tz] = waves_fdm_1d(def,sigma,plot_flag,order,true);
    errors(i) = max(e);
    errors_tz(i) = max(e_tz);
end
hp = h.^order;
figure(2)
loglog(h,errors,'o',h,errors_tz,'o',h,hp);
xlabel("h");
ylabel("|e|_{\infty}");
title("Truncation Error");
str = sprintf("h^{%d}",order);
legend("error","error_{tz}",str);
