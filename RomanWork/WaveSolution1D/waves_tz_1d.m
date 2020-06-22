% waves_tz_1d.m
% twilight zone method for code verification

function u = waves_tz_1d(def,sigma,plot_flag)
    dx = (def.b-def.a)/def.N;
    nt = def.t_f/(sigma*dx/def.c);
    nt = ceil(nt);
    dt = def.t_f/nt;
    sigma = def.c*dt/dx;
    
    % v(x,t): our "exact soln"
    function val = v(x,t)
        val = sin(x)*sin(t);
    end

    % v_t(x,t)
    function val = v_t(x,t)
        val = cos(t)*sin(x);
    end

    % v_x(x,t)
    function val = v_x(x,t)
        val = cos(x)*sin(t);
    end

    % v_xx(x,t)
    function val = v_xx(x,t)
        val = -sin(x)*sin(t);
    end

    % h = (u_e)_tt-c^2(u_e)_xx
    function val = h(x,t)
        val = (def.c^2-1)*sin(x)*sin(t);
    end
    
    x = linspace(def.a,def.b,def.N+1);
    
    % initial conditions
    unm1 = zeros(1,def.N+1);
    for j = 1:def.N+1
        unm1(j) = v(x(j),0);
    end
    % first time step
    un = zeros(1,def.N+1);
    un(1) = def.left(dt);
    for j = 2:def.N
        un(j) =  v(x(j),0)+dt*v_t(x(j),0)+(dt)^2/2*(def.c^2*v_xx(x(j),0)+h(x(j),0));
%         un(j) = v(x(j),dt);
    end
    un(def.N+1) = (4*un(def.N)-un(def.N-1)-2*def.right(dt)*dx)/3;
    % time loop - for n >= 2
    n = 2;
    while n*dt <= def.t_f
        unp1 = zeros(1,def.N+1);
        unp1(1) = def.left(n*dt);
%         unp1(1) = v(0,n*dt);
        for j = 2:def.N
            unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1))+dt^2*h(x(j),n*dt);
        end
        unp1(def.N+1) = (4*unp1(def.N)-unp1(def.N-1)-2*def.right(n*dt)*dx)/3;
%         unp1(def.N+1) = (4*unp1(def.N)-unp1(def.N-1)-2*v_x(x(def.N+1),n*dt)*dx)/3;
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
end

