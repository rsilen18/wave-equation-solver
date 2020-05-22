% driver code for waves_fdm_1d.m
N = 100;    % resolution in x-dimension
c = 1;      % wave speed
t_f = 10;   % final time
sigma = 1.0005; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
plot_flag = true;    % flag for turning on plot animation

% Original problem
[u,e] = waves_fdm_1d(N,sigma,c,t_f,icase,plot_flag);

% Convergence study
