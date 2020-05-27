% driver code for waves_fdm_1d.m
N = 100;    % resolution in x-dimension
c = 1;      % wave speed
t_f = 10;   % final time
sigma = .5; % CFL condition: sigma <= 1
icase = 1;  % flag for problem definition
plot_flag = true;    % flag for turning on plot animation

% Original problem
[u,e] = waves_fdm_1d(N,sigma,c,t_f,icase,plot_flag);

% % Convergence study
% icase = 5;
% plot_flag = false;
% dx = linspace(1e-4,1e-1,10);
% N = round(ones(1,10)./dx);
% errors = zeros(1,10);
% for i = 1:10
%     [u,errors(i)] = waves_fdm_1d(N(i),sigma,c,t_f,icase,plot_flag);
% end
% loglog(dx,errors);
% xlabel("dx");
% ylabel("|e|_{\infty}");
% title("error vs. dx");