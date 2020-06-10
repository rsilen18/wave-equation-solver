% 1-D Finite Difference Method Wave Eqn Solver
% BCs: u(a,t) = l(t), u_x(b,t) = r(t)
% ICs: u(x,0) = u_0(x), u_t(x,0) = v_0(x)

function [u,e] = waves_fdm_1d(def,sigma,plot_flag,order)
    dx = (def.b-def.a)/def.N;
    nt = def.t_f/(sigma*dx/def.c);
    nt = ceil(nt);
    dt = def.t_f/nt;
    sigma = def.c*dt/dx;
    
    % to account for index offset in 4th order
    function index = j_(i)
        index = i+3;
    end
    
    if (order == 2)
        x = linspace(def.a,def.b,def.N+1);
        % initial conditions
        unm1 = zeros(1,def.N+1);
        for j = 1:def.N+1
            unm1(j) = def.f(x(j));
        end
        unm1(def.N+1) = (4*unm1(def.N)-unm1(def.N-1)-2*def.right(0)*dx)/3;
        un = zeros(1,def.N+1);
        un(1) = def.left(dt);
        for j = 2:def.N
           un(j) =  (1-sigma^2)*unm1(j)+dt*def.g(x(j))+sigma^2/2*(unm1(j-1)+unm1(j+1));
        end
        un(def.N+1) = (4*un(def.N)-un(def.N-1)-2*def.right(dt)*dx)/3;
        % time loop
        n = 2;
        while n*dt <= def.t_f
            unp1 = zeros(1,def.N+1);
            unp1(1) = def.left(n*dt);
            for j = 2:def.N
                unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1));
            end
            unp1(def.N+1) = (4*unp1(def.N)-unp1(def.N-1)-2*def.right(n*dt)*dx)/3;
            % optionally plot
            if (plot_flag)
               plot(x,unp1);
               xlim([def.a def.b]);
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
        e = abs(u-waves_analytic_1d(def.f,def.g,def.c,x,def.t_f));
        
    elseif (order == 4)
        x = linspace(def.a-2*dx,def.b+2*dx,def.N+5);
        
        % first time step
        unm1 = zeros(def.N+5);
        un = zeros(def.N+5);
        for i = 1:def.N+5
            unm1(i) = waves_analytic_1d(def.f,def.g,def.c,x(i),0);
            un(i) = waves_analytic_1d(def.f,def.g,def.c,x(i),dt);
        end
        % remaining steps
        n = 2;
        while n*dt <= def.t_f
            unp1 = zeros(1,def.N+5);
            unp1(j_(-2)) = waves_analytic_1d(def.f,def.g,def.c,x(j_(-2)),n*dt);
            unp1(j_(-1)) = waves_analytic_1d(def.f,def.g,def.c,x(j_(-1)),n*dt);
            unp1(j_(def.N+1)) = waves_analytic_1d(def.f,def.g,def.c,x(j_(def.N+1)),n*dt);
            unp1(j_(def.N+2)) = waves_analytic_1d(def.f,def.g,def.c,x(j_(def.N+2)),n*dt);
            for j = 0:def.N
                order2 = un(j_(j+1))-2*un(j_(j))+un(j_(j-1));
                order4 = un(j_(j+2))-4*un(j_(j+1))+6*un(j_(j))-4*un(j_(j-1))+un(j_(j-2));
                unp1(j_(j)) = 2*un(j_(j))-unm1(j_(j))+sigma^2*(order2-1/12*order4)+1/12*sigma^4*order4;
            end
            % optionally plot
            if (plot_flag)
               plot(x,unp1);
               xlim([def.a def.b]);
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
        u = unp1(j_(0):j_(def.N));
        e = abs(u-waves_analytic_1d(def.f,def.g,def.c,x(j_(0):j_(def.N)),def.t_f));
    end
end
