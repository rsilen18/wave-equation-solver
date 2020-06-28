% waves_tz_1d.m
% twilight zone method for code verification

function [u,e] = waves_tz_1d(def,sigma,plot_flag,order)
    dx = (def.b-def.a)/def.N;
    nt = def.t_f/(sigma*dx/def.c);
    nt = ceil(nt);
    dt = def.t_f/nt;
    sigma = def.c*dt/dx;
    
    x = linspace(def.a-dx*order/2, def.b+dx*order/2, def.N+order+1);
    unm1 = ICs(def,x);
    un = first_time_step(def,sigma,x,order,dt,unm1);
    un = BCs(def,sigma,x,order,dx,dt,un,1);
    n = 2;
    while n*dt <= def.t_f
        unp1 = main_time_step(def,sigma,x,order,dt,unm1,un,n);
        unp1 = BCs(def,sigma,x,order,dx,dt,unp1,n);
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
    u = unp1(1+order/2:size(x,2)-order/2);
    exact = v(x(1+order/2:size(x,2)-order/2),def.t_f);
    e = abs(u-exact);
    
end

%%% forcing functions

% v(x,t): our "exact soln"
function val = v(x,t)
    val = sin(x).*sin(t);
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

function val = v_xxt(x,t)
    val = -sin(x)*cos(t);
end

function val = v_4x(x,t)
    val = sin(x)*sin(t);
end

% v_xxxxt(x,t)
function val = v_4xt(x,t)
    val = sin(x)*cos(t);
end

function val = v_6x(x,t)
    val = -sin(x)*sin(t);
end

% h = v_tt-c^2*v_xx
function val = h(def,x,t)
    val = (def.c^2-1)*sin(x)*sin(t);
end

function val = h_t(def,x,t)
    val = (def.c^2-1)*sin(x)*cos(t);
end

function val = h_xx(def,x,t)
    val = (1-def.c^2)*sin(x)*sin(t);
end

function val = h_tt(def,x,t)
    val = (1-def.c^2)*sin(x)*sin(t);
end

function val = h_xxt(def,x,t)
    val = (1-def.c^2)*sin(x)*cos(t);
end

function val = h_3t(def,x,t)
    val = (1-def.c^2)*sin(x)*cos(t);
end

function val = h_4x(def,x,t)
    val = (def.c^2-1)*sin(x)*sin(t);
end

function val = h_4t(def,x,t)
    val = (def.c^2-1)*sin(x)*sin(t);
end

function val = h_xxtt(def,x,t)
    val = (def.c^2-1)*sin(x)*sin(t);
end

%%% set up functions

% initial conditions
function unm1 = ICs(def,x)
    n = size(x,2);
    unm1 = zeros(1,n);
    for j = 1:n
        unm1(j) = v(x(j),0);
    end
end

% fill in boundary conditions
% code for Dirichlet left BC, Neumann right BC
% n = timestep
function un = BCs(def,sigma,x,order,dx,dt,un,n)
    m = size(x,2);
    ja = 1+order/2;
    jb = m-order/2;
    if (order == 2)
        un(ja-1) = 2*un(ja)-un(ja+1);
        un(jb+1) = un(jb-1)+2*dx*v_x(x(jb),n*dt);
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
        A = [f1-f0, f2-f0];
        b = -1*f0;
        u = A\b;
        un(ja-1) = u(1);
        un(ja-2) = u(2);

        % Neumann right BC - "way 1"
        A = [2/3 -1/12; -2 1];
        b = [2/3*un(jb-1)-1/12*un(jb-2)+dx*v_x(x(jb),n*dt);
            -2*un(jb-1)+un(jb-2)+2*dx/sigma^2*(v_x(x(jb),(n+1)*dt)-...
            2*v_x(x(jb),n*dt)+v_x(x(jb),(n-1)*dt))];
        u = A\b;
        un(jb+1) = u(1);
        un(jb+2) = u(2);
    elseif (order == 6)
        % Dirichlet left BC - discrete delta fcn
        % u(1) = un(ja-1)
        % u(2) = un(ja-2)
        % u(3) = un(ja-3)
%             f = @(u) [def.c^2/dx^2*(1/90*un(ja+3)-3/20*un(ja+2)+3/2*un(ja+1)-49/18*un(ja)+3/2*u(1)-3/20*u(2)+1/90*u(3))-def.left(n*dt);
%                 def.c^4/dx^4*(-1/6*un(ja+3)+2*un(ja+2)-13/2*un(ja+1)+28/3*un(ja)-13/2*u(1)+2*u(2)-1/6*u(3));
%                 def.c^6/dx^6*(un(ja+3)-6*un(ja+2)+15*un(ja+1)-20*un(ja)+15*u(1)-6*u(2)+u(3))];
        f = @(u) [(1/90*un(ja+3)-3/20*un(ja+2)+3/2*un(ja+1)-49/18*un(ja)+3/2*u(1)-3/20*u(2)+1/90*u(3)-v(x(ja),n*dt));
            (-1/6*un(ja+3)+2*un(ja+2)-13/2*un(ja+1)+28/3*un(ja)-13/2*u(1)+2*u(2)-1/6*u(3));
            (un(ja+3)-6*un(ja+2)+15*un(ja+1)-20*un(ja)+15*u(1)-6*u(2)+u(3))];
        f0 = f([0;0;0]);
        f1 = f([1;0;0]);
        f2 = f([0;1;0]);
        f3 = f([0;0;1]);
        A = [f1-f0, f2-f0, f3-f0];
        b = -1*f0;
        u = A\b;
        un(ja-1) = u(1);
        un(ja-2) = u(2);
        un(ja-3) = u(3);

        % Neumann right BC - discrete delta fcn
        d = 30; % unknown denominator of coefficient(?)
%             f = @(u) [1/dx*(u(1)-un(jb-1)-1/6*(u(2)-2*u(1)+2*un(jb-1)-un(jb-2))+...
%                 1/d*(u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3)))-def.right(n*dt);
%                 def.c^2/dx^3*(u(2)-2*u(1)+2*un(jb-1)-un(jb-2)+...
%                 1/d*(u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3)));
%                 def.c^4/dx^5*(u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3))];
        f = @(u) [(u(1)-un(jb-1)-1/6*(u(2)-2*u(1)+2*un(jb-1)-un(jb-2))+...
            1/d*(u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3)))-v_x(x(jb),n*dt);
            (u(2)-2*u(1)+2*un(jb-1)-un(jb-2)+...
            1/d*(u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3)));
            (u(3)-4*u(2)+5*u(1)-5*un(jb-1)+4*un(jb-2)-un(jb-3))];
        f0 = f([0;0;0]);
        f1 = f([1;0;0]);
        f2 = f([0;1;0]);
        f3 = f([0;0;1]);
        A = [f1-f0, f2-f0, f3-f0];
        b = -1*f0;
        u = A\b;
        un(jb+1) = u(1);
        un(jb+2) = u(2);
        un(jb+3) = u(3);
    end
end

% first time step
function un = first_time_step(def,sigma,x,order,dt,unm1)
    n = size(x,2);
    un = zeros(1,n);
    for j = 1+order/2:n-order/2
        if (order >= 2)
            un(j) =  unm1(j)+dt*v_t(x(j),0)+dt^2/2*(def.c^2*v_xx(x(j),0)+h(def,x(j),0));
        end
        if (order >= 4)
            un(j) = un(j) + ...
                    dt^3/6*(def.c^2*v_xxt(x(j),0)+h_t(def,x(j),0)) + ...
                    dt^4/24*(def.c^4*v_4x(x(j),0)+def.c^4*h_xx(def,x(j),0)+h_tt(def,x(j),0));
        end
        if (order >= 6)
            un(j) = un(j) + ...
                    dt^5/120*(def.c^4*v_4xt(x(j),0)+def.c^2*h_xxt(def,x(j),0)+h_3t(def,x(j),0)) + ...
                    dt^6/720*(def.c^6*v_6x(x(j),0)+def.c^4*h_4x(def,x(j),0)+h_xxtt(def,x(j),0));
        end
    end
end

% time steps for n >= 2
function unp1 = main_time_step(def,sigma,x,order,dt,unm1,un,n)
    m = size(x,2);
    unp1 = zeros(1,m);
    for j = 1+order/2:m-order/2
        unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1))+h(def,x(j),n*dt);
        if (order >= 4)
            unp1(j) = unp1(j) + (sigma^4-sigma^2)/12*(un(j+2)-4*un(j+1)+6*un(j)-4*un(j-1)+un(j-2)) + ...
                def.c^2*dt^4/12*h_xx(def,x(j),n*dt) + dt^4/12*h_tt(def,x(j),n*dt);
        end
        if (order >= 6)
%                 unp1(j) = unp1(j) + (4*sigma^2+sigma^6)/360*...
%                 (un(j+3)-6*un(j+2)+15*un(j+1)-20*un(j)+15*un(j-1)-6*un(j-2)+un(j-3));
            unp1(j) = unp1(j) + (sigma^2/90-sigma^4/72+sigma^6/720)*...
                (un(j+3)-6*un(j+2)+15*un(j+1)-20*un(j)+15*un(j-1)-6*un(j-2)+un(j-3)) + ...
                def.c^2*dt^6/360*(h_xx(def,x(j),n*dt)+h_xxtt(def,x(j),0)+h_4t(def,x(j),0));
        end
    end
end
