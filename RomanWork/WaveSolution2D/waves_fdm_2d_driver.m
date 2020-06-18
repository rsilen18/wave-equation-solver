% driver code for waves_fdm_2d.m
sigma = .5; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
plot_flag = true;    % flag for turning on plot animation

% problem definition
[def.a_x,def.b_x,def.a_y,def.b_y,def.c,def.N,def.t_f,...
    def.f,def.g,def.left,def.right,def.bottom,def.top]...
    = waves_fdm_2d_defs(icase);

% Original problem
[u,e] = waves_fdm_2d(def,sigma,plot_flag);
