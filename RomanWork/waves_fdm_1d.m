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
    
    % to account for index offset in 4th order
    function index = j_(i)
        index = i+3;
    end
    
    if (order == 2)
        x = linspace(a,b,N+1);
        % initial conditions
        unm1 = zeros(1,N+1);
        for j = 1:N+1
            unm1(j) = u_0(x(j));
        end
        unm1(N+1) = (4*unm1(N)-unm1(N-1)-2*right(0)*dx)/3;
        un = zeros(1,N+1);
        un(1) = left(dt);
        for j = 2:N
           un(j) =  (1-sigma^2)*unm1(j)+dt*v_0(x(j))+sigma^2/2*(unm1(j-1)+unm1(j+1));
        end
        un(N+1) = (4*un(N)-un(N-1)-2*right(dt)*dx)/3;
        % time loop
        n = 2;
        while n*dt <= t_f
            unp1 = zeros(1,N+1);
            unp1(1) = left(n*dt);
            for j = 2:N
                unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1));
            end
            unp1(N+1) = (4*unp1(N)-unp1(N-1)-2*right(n*dt)*dx)/3;
            % optionally plot
            if (plot_flag)
               plot(x,unp1);
               xlim([a b]);
               ylim([-2 2]);
               xlabel("x step");
               ylabel("displacement: u(x,t)");
               str = sprintf("1-D Wave at t=%.2f",n*dt);
               title(str);
               pause(.0001);
            end
            unm1 = un;
            un = unp1;
            n = n+1;
        end
        % assignment of u at final time step
        u = unp1;
        e = abs(u-waves_analytic_1d(u_0,v_0,c,x,t_f));
        
    elseif (order == 4)
        x = linspace(a-2*dx,b+2*dx,N+5);
        
        % first time step
        unm1 = waves_analytic_1d(u_0,v_0,c,x,0);
        un = waves_analytic_1d(u_0,v_0,c,x,dt);
        % remaining steps
        n = 2;
        while n*dt <= t_f
            unp1 = zeros(1,N+5);
            unp1(j_(-2)) = waves_analytic_1d(u_0,v_0,c,x(j_(-2)),n*dt);
            unp1(j_(-1)) = waves_analytic_1d(u_0,v_0,c,x(j_(-1)),n*dt);
            unp1(j_(N+1)) = waves_analytic_1d(u_0,v_0,c,x(j_(N+1)),n*dt);
            unp1(j_(N+2)) = waves_analytic_1d(u_0,v_0,c,x(j_(N+2)),n*dt);
            for j = 0:N
                order2 = un(j_(j+1))-2*un(j_(j))+un(j_(j-1));
                order4 = un(j_(j+2))-4*un(j_(j+1))+6*un(j_(j))-4*un(j_(j-1))+un(j_(j-2));
                unp1(j_(j)) = 2*un(j_(j))-unm1(j_(j))+sigma^2*(order2-1/12*order4)+1/12*sigma^4*order4;
            end
            % optionally plot
            if (plot_flag)
               plot(x,unp1);
               xlim([a b]);
               ylim([-2 2]);
               xlabel("x step");
               ylabel("displacement: u(x,t)");
               str = sprintf("1-D Wave at t=%.2f",n*dt);
               title(str);
               pause(.0001);
            end
            unm1 = un;
            un = unp1;
            n = n+1;
        end
        % assignment of u at final time step
        u = unp1(j_(0):j_(N));
        e = abs(u-waves_analytic_1d(u_0,v_0,c,x(j_(0):j_(N)),t_f));
    end
end
