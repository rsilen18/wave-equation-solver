% driver code for waves_fdm_2d.m
sigma_x = .5; % CFL condition: sigma <= 1
sigma_y = .5; 
icase = 1;  % flag for problem definition
order = 2;
plot_flag = true;    % flag for turning on plot animation

% problem definition
[def.a_x,def.b_x,def.a_y,def.b_y,def.c,def.N,def.t_f,...
    def.f,def.g,def.left,def.right,def.bottom,def.top]...
    = waves_fdm_2d_defs(icase);

% Original problem
% u = waves_fdm_2d(def,sigma_x,sigma_y,plot_flag,order);

% Convergence study
icase = 1;
plot_flag = false;
h = [.1 .01 ];
N = [10 100 ];
errors_tz = zeros(1,2);
for i = 1:2
    disp(i);
    def.N = N(i);
    [u_tz,e_tz] = waves_tz_2d(def,sigma_x,sigma_y,plot_flag,order);
    errors_tz(i) = max(max(e_tz));
end
h2 = h.^2;
h3 = h.^3;
h4 = h.^4;
h5 = h.^5;
h6 = h.^6;
figure(2)
loglog(h,errors_tz,'o',h,h2,h,h3,h,h4,h,h5,h,h6);
xlabel("h");
ylabel("|e|_{\infty}");
title("Truncation Error");
legend("error_{tz}","h^2","h^3","h^4","h^5","h^6");