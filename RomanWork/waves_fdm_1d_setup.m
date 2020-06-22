% waves_fdm_1d_setup.m
% set up functions for waves_fdm_1d.m

% initial conditions
    function unm1 = ICs(x)
        n = size(x,2);
        unm1 = zeros(1,n);
        for j = 1:n
            unm1(j) = def.f(x(j));
        end
    end

    % fill in boundary conditions
    % code for Dirichlet left BC, Neumann right BC
    % n = timestep
    function un = BCs(x,un,n)
        m = size(x,2);
        ja = 1+order/2;
        jb = m-order/2;
        if (order == 2)
            un(ja-1) = 2*un(ja)-un(ja+1);
            un(jb+1) = un(jb-1);
        elseif (order == 4)
            % Dirichlet left BC - discrete delta fcn 
            % u(1) = un(ja-1)
            % u(2) = un(ja-2)
            f = @(u) [def.c^2/dx^2*(un(ja+1)-2*un(ja)+u(1) - ...
                1/12*(un(ja+2)-4*un(ja+1)+6*un(ja)-4*u(1)+u(2)));
                def.c^4/dx^4*(un(ja+2)-4*un(ja+1)+6*un(ja)-4*u(1)+u(2))];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0; f2-f0];
            b = -1*f0;
            u = A\b;
            un(ja-1) = u(1);
            un(ja-2) = u(2);
            
            % Neumann right BC - "way 1"
            A = [2/3 -1/12; -2 1];
            b = [2/3*un(jb-1)-un(jb-2)+dx*def.right(n*dt);
                -2*un(jb-1)+un(jb-2)+2*dx^3/(def.c^2*dt^2)*(def.right((n+1)*dt)-...
                2*def.right(n*dt)+def.right((n-1)*dt))];
            u = A\b;
            un(jb+1) = u(1);
            un(jb+2) = u(2);
        end
    end

    % first time step
    function un = first_time_step(x)
        n = size(x,2);
        un = zeros(1,n);
        for j = 1+order/2:n-order/2
            if (order == 2)
                un(j) =  (1-sigma^2)*unm1(j)+dt*def.g(x(j))+sigma^2/2*(unm1(j-1)+unm1(j+1));
            elseif (order == 4)
                un(j) = unm1(j) + ...
                        dt*def.g(x(j)) + ...
                        sigma^2/2*(unm1(j+1)-2*unm1(j)+unm1(j-1)-1/12*(unm1(j+2)-4*unm1(j+1)+6*unm1(j)-4*unm1(j-1)+unm1(j-2))) + ...
                        dt*sigma^2/6*(def.g(x(j+1))-2*def.g(x(j))+def.g(x(j-1))-1/12*(def.g(x(j+2))-4*def.g(x(j+1))+6*def.g(x(j))-4*def.g(x(j-1))+def.g(x(j-2)))) + ...
                        sigma^4/24*(unm1(j+2)-4*unm1(j+1)+6*unm1(j)-4*unm1(j-1)+unm1(j-2));
            end
        end
    end

    % time steps for n >= 2
    function unp1 = main_time_step(x)
        n = size(x,2);
        unp1 = zeros(1,n);
        for j = 1+order/2:n-order/2
            unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1));
            if (order == 4)
                unp1(j) = unp1(j) + (sigma^4-sigma^2)/12*(un(j+2)-4*un(j+1)+6*un(j)-4*un(j-1)+un(j-2));
            end
        end
    end