% 1-D Finite Difference Method Wave Eqn Solver
% BCs: u(a,t) = l(t), u_x(b,t) = r(t)
% ICs: u(x,0) = u_0(x), u_t(x,0) = v_0(x)

function [u,e] = waves_fdm_1d(N,sigma,c,t_f,icase,plot_flag,order)
    [a,b,u_0,v_0,left,right] = waves_fdm_1d_defs(icase);
    dx = (b-a)/N;
    nt = t_f/(sigma*dx/c);
    nt = ceil(nt);
    dt = t_f/nt;
    sigma = c*dt/dx;
    
    if (order == 2)
        x = linspace(a,b,N+1);
        % initial conditions
        unm1 = zeros(1,N+1);
        for i = 1:N+1
            unm1(i) = u_0(x(i));
        end
        unm1(N+1) = (4*unm1(N)-unm1(N-1)-2*right(0)*dx)/3;
        un = zeros(1,N+1);
        un(1) = left(dt);
        for i = 2:N
           un(i) =  (1-sigma^2)*unm1(i)+dt*v_0(x(i))+sigma^2/2*(unm1(i-1)+unm1(i+1));
        end
        un(N+1) = (4*un(N)-un(N-1)-2*right(dt)*dx)/3;
        % time loop
        j = 2;
        while j*dt <= t_f
            unp1 = zeros(1,N+1);
            unp1(1) = left(j*dt);
            for i = 2:N
                unp1(i) = 2*un(i)-unm1(i)+sigma^2*(un(i+1)-2*un(i)+un(i-1));
            end
            unp1(N+1) = (4*unp1(N)-unp1(N-1)-2*right(j*dt)*dx)/3;
            % optionally plot
            if (plot_flag)
               plot(x,unp1);
               xlim([a b]);
               ylim([-2 2]);
               xlabel("x step");
               ylabel("displacement: u(x,t)");
               str = sprintf("1-D Wave at t=%.2f",j*dt);
               title(str);
               pause(.0001);
            end
            unm1 = un;
            un = unp1;
            j = j+1;
        end
        % assignment of u at final time step
        u = unp1;
        e = abs(u-waves_analytic_1d(u_0,v_0,c,x,t_f));
        
    elseif (order == 4)
        x = linspace(a-2*dx,b+2*dx,N+4);
        j = 1;
        % first time step
        u1 = 2;
    end
end
