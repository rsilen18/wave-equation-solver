% 2-D Finite Difference Method Wave Eqn Solver
% BCs: u(0,y,t)=u(a,y,t) = 0        0 <= y <= b
%      u(x,0,t)=u(x,b,t) = 0        0 <= x <= a
% ICs: u(x,y,0) = f(x,y), u_t(x,y,0) = g(x,y)

function [u,e] = waves_fdm_2d(def,sigma,plot_flag)
    % we let dx = dy
    dx = (def.b-def.a)/def.N;
    nt = def.t_f/(sigma*dx/def.c);
    nt = ceil(nt);
    dt = def.t_f/nt;
    sigma = def.c*dt/dx;
    
    function unp1 = fill_BCs(def,x,y,t)
        unp1 = zeros(def.N,def.N);
        for i_ = 1:def.N
            unp1(i_,1) = def.bottom(x(i_),t);
            unp1(i_,def.N) = def.top(x(i_),t);
        end
        for j_ = 1:def.N
            unp1(1,j_) = def.left(y(j_),t);
            unp1(def.N,j_) = def.right(y(j_),t);
        end
    end

    x = linspace(def.a_x,def.b_x,def.N);
    y = linspace(def.a_y,def.b_y,def.N);
    % initial conditions
    unm1 = zeros(def.N,def.N);
    for i = 1:def.N
        for j = 1:def.N
            unm1(i,j) = def.f(x(i),y(j));
        end
    end
    % first time step
    un = fill_BCs(def,x,y,dt);
    for i = 2:def.N-1
        for j = 2:def.N-1
            un(i,j) = (1-sigma^2)*unm1(i,j)+dt*def.g(x(i),y(j))...
                +sigma^2/4*(unm1(i-1,j)+unm1(i+1,j)+unm1(i,j-1)...
                +unm1(i,j+1));
        end
    end
    % time loop
    n = 2;
    while n*dt <= def.t_f
        unp1 = fill_BCs(def,x,y,n*dt);
        for i = 2:def.N-1
            for j = 2:def.N-1
                unp1(i,j) = 2*un(i,j)-unm1(i,j)+sigma^2*(un(i-1,j)...
                    +un(i+1,j)+un(i,j-1)+un(i,j+1)-4*un(i,j));
            end
        end
        % optionally plot
        if (plot_flag)
           surf(x,y,unp1);
           xlim([def.a_x def.b_x]);
           ylim([def.a_y def.b_y]);
           zlim([-1 1]);
           xlabel("x step");
           ylabel("y step");
           zlabel("displacement u(x,y,t)");
           str = sprintf("2-D Wave at t=%.2f",n*dt);
           title(str);
           pause(.0001);
        end
        unm1 = un;
        un = unp1;
        n = n+1;
    end
    % assignment of u at final time step
    u = unp1;
    e = 0;
%     e = abs(u-waves_analytic_1d(def.f,def.g,def.c,x,def.t_f));
end
