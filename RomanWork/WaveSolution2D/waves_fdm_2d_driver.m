% driver code for waves_fdm_2d.m
cfl = .5;   % CFL condition: cfl in (0,1)
icase = 1;  % flag for problem definition
order = 4;
plot_flag = false;    % flag for turning on plot animation

% problem definition
[def.a_x,def.b_x,def.a_y,def.b_y,def.c,def.N,def.t_f,...
    def.f,def.g,def.left,def.right,def.bottom,def.top]...
    = waves_fdm_2d_defs(icase);

% Original problem
% u = waves_fdm_2d(def,cfl,plot_flag,order);

% [u_tz,e_tz] = waves_tz_2d(def,cfl,plot_flag,order);

% Convergence study
icase = 1;
plot_flag = false;
h = [.1 .05 .04 .025 .02 .125 .01];
N = 1./h;
errors_tz = zeros(1,size(h,2));
for i = 1:size(h,2)
    disp(i);
    def.N = N(i);
    [u_tz,e_tz] = waves_tz_2d(def,cfl,plot_flag,order);
    errors_tz(i) = max(max(e_tz));
%     figure(i+1)
    % plot error surface
%     x = linspace(def.a_x,def.b_x,def.N+1);
%     y = linspace(def.a_y,def.b_y,def.N+1);
%     surf(x,y,e_tz);
%     xlabel("x");
%     ylabel("y");
%     str = sprintf("Error surface plot for h=%f",h(i));
%     title(str);
end
hp = h.^order;
figure(size(h,2)+2)
loglog(h,errors_tz,'o',h,hp);
xlabel("h");
ylabel("|e|_{\infty}");
title("Truncation Error");
str = sprintf("h^{%d}",order);
legend("error_{tz}",str);