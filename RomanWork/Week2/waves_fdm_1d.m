% 1-D Finite Difference Method Wave Eqn Solver
% BCs: u(a,t) = l(t), u_x(b,t) = 0
% ICs: u(x,0) = u_0(x), u_t(x,0) = v_0(x)

function [u,e] = waves_fdm_1d(N,sigma,c,t_f,icase,plot_flag)
    [a,b,u_0,v_0,left] = waves_fdm_1d_defs(icase);
    dx = (b-a)/N;
    nt = t_f/(sigma*dx/c);
    nt = ceil(nt);
    dt = t_f/nt;
    x = linspace(a,b,N+1);
    
    % initial conditions
    unm1 = zeros(1,N+1);
    for i = 1:N+1
        unm1(i) = u_0(x(i));
    end
    ghost = unm1(N);
    un = zeros(1,N+1);
    un(1) = left(dt);
    for i = 2:N
       un(i) =  (1-sigma^2)*unm1(i)+dt*v_0(x(i))+sigma^2/2*(unm1(i-1)+unm1(i+1));
    end
    un(N+1) = (1-sigma^2)*unm1(i)+dt*v_0(x(i))+sigma^2/2*(unm1(i-1)+ghost);
    ghost = un(N);
    % time loop
    j = 2;
    while j*dt <= t_f
        unp1 = zeros(1,N+1);
        unp1(1) = left(j*dt);
        for i = 2:N
            unp1(i) = 2*un(i)-unm1(i)+sigma^2*(un(i+1)-2*un(i)+un(i-1));
        end
        unp1(N+1) = 2*un(i)-unm1(i)+sigma^2*(ghost-2*un(i)+un(i-1));
        ghost = unp1(N);
        % optionally plot
        if (plot_flag)
           plot(x,unp1);
           xlim([a b]);
           ylim([-2 2]);
           xlabel("x step");
           ylabel("displacement: u(x,t)");
           str = sprintf("1-D Wave at t=%.2f",j*dt);
           title(str);
           pause(.001);
        end
        unm1 = un;
        un = unp1;
        j = j+1;
    end
    % assignment of u at final time step
    u = unp1;
    e = max(abs(waves_analytic_1d(u_0,v_0,a,b,c,x,t_f)-u));
end
