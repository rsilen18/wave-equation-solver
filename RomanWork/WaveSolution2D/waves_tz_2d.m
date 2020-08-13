% waves_tz_2d.m
% 2-D twilight zone verification
% BCs: u(0,y,t)=u(a,y,t) = 0        0 <= y <= b
%      u(x,0,t)=u(x,b,t) = 0        0 <= x <= a
% ICs: u(x,y,0) = f(x,y), u_t(x,y,0) = g(x,y)

function [u,e] = waves_tz_2d(def,cfl,plot_flag,order)
    dx = (def.b_x-def.a_x)/def.N;
    dy = (def.b_y-def.a_y)/def.N;
    dt = cfl*dx*dy/def.c*sqrt(1/(dx^2+dy^2));
    nt = ceil(def.t_f/dt);
    dt = def.t_f/nt;
    sigma_x = def.c*dt/dx;
    sigma_y = def.c*dt/dy;

    x = linspace(def.a_x-dx*order/2,def.b_x+dx*order/2,def.N+order+1);
    y = linspace(def.a_y-dy*order/2,def.b_y+dy*order/2,def.N+order+1);
    
    unm1 = ICs(def,x,y);
    un = first_time_step(def,sigma_x,sigma_y,x,y,order,dt,unm1);
    un = BCs(def,x,y,order,dx,dy,un,1,dt);
    n = 2;
    while n*dt <= def.t_f
        unp1 = main_time_step(def,sigma_x,sigma_y,x,y,order,dt,dx,dy,unm1,un,n-1);
        unp1 = BCs(def,x,y,order,dx,dy,unp1,n,dt);
        % optionally plot
        if (plot_flag)
%            surf(x,y,unp1);
%            shading interp;
%            xlim([def.a_x def.b_x]);
%            ylim([def.a_y def.b_y]);
%            zlim([-3 3]);
%            xlabel("x step");
%            ylabel("y step");
%            zlabel("displacement u(x,y,t)");
%            str = sprintf("2-D Wave at t=%.2f",n*dt);
%            title(str);
%            pause(.0001);
%             if (mod(n,20) == 0)
                u = unp1(1+order/2:size(x,2)-order/2,1+order/2:size(y,2)-order/2);
                exact = zeros(def.N+order+1,def.N+order+1);
                for i = 1:def.N+order+1
                    for j = 1:def.N+order+1
                        exact(i,j) = v(x(i),y(j),n*dt);
                    end
                end
                e = abs(u-exact(1+order/2:size(x,2)-order/2,1+order/2:size(y,2)-order/2));
                surf(x(1+order/2:size(x,2)-order/2),y(1+order/2:size(y,2)-order/2),e);
                shading interp;
                str = sprintf("Error at t=%.4f",n*dt);
                xlabel("x step");
                ylabel("y step");
                title(str);
                drawnow;
                pause;
%             end
        end
        unm1 = un;
        un = unp1;
        n = n+1;
    end
    % assignment of u at final time step
    u = unp1(1+order/2:size(x,2)-order/2,1+order/2:size(y,2)-order/2);
    exact = zeros(def.N+order+1,def.N+order+1);
    for i = 1:def.N+order+1
        for j = 1:def.N+order+1
            exact(i,j) = v(x(i),y(j),def.t_f);
        end
    end
    e = abs(u-exact(1+order/2:size(x,2)-order/2,1+order/2:size(y,2)-order/2));
end

%%% forcing functions

% v(x,y,t): our "exact soln"
function val = v(x,y,t)
    val = sin(x)*sin(y)*sin(t);
%     val = sin(t);
end

% v_t(x,y,t)
function val = v_t(x,y,t)
    val = sin(x)*sin(y)*cos(t);
%     val = cos(t);
end

% v_x(x,y,t)
function val = v_x(x,y,t)
    val = cos(x)*sin(y)*sin(t);
%     val = 0;
end

% v_y(x,y,t)
function val = v_y(x,y,t)
    val = sin(x)*cos(y)*sin(t);
%     val = 0;
end

function val = v_xx(x,y,t)
    val = -sin(x)*sin(y)*sin(t);
%     val = 0;
end

function val = v_yy(x,y,t)
    val = -sin(x)*sin(y)*sin(t);
%     val = 0;
end

function val = v_tt(x,y,t)
    val = -sin(x)*sin(y)*sin(t);
end

function val = v_xxt(x,y,t)
    val = -sin(x)*sin(y)*cos(t);
end

function val = v_yyt(x,y,t)
    val = -sin(x)*sin(y)*cos(t);
end

function val = v_4x(x,y,t)
    val = sin(x)*sin(y)*sin(t);
end

function val = v_xxyy(x,y,t)
    val = sin(x)*sin(y)*sin(t);
end

function val = v_4y(x,y,t)
    val = sin(x)*sin(y)*sin(t);
end

% h = v_tt-c^2*(v_xx+v_yy)
function val = h(def,x,y,t)
%     val = (2*def.c^2-1)*sin(x)*sin(y)*sin(t);
    val = v_tt(x,y,t)-def.c^2*(v_xx(x,y,t)+v_yy(x,y,t));
end

function val = h_t(def,x,y,t)
    val = (2*def.c^2-1)*sin(x)*sin(y)*cos(t);
%     val = -cos(t);
end

function val = h_xx(def,x,y,t)
    val = (1-2*def.c^2)*sin(x)*sin(y)*sin(t);
end

function val = h_yy(def,x,y,t)
    val = (1-2*def.c^2)*sin(x)*sin(y)*sin(t);
end

function val = h_tt(def,x,y,t)
    val = (1-2*def.c^2)*sin(x)*sin(y)*sin(t);
end

%%% set up functions

% initial conditions
function unm1 = ICs(def,x,y)
    nx = size(x,2);
    ny = size(y,2);
    unm1 = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            unm1(i,j) = v(x(i),y(j),0);
        end
    end
end

% first time step
function un = first_time_step(def,sigma_x,sigma_y,x,y,order,dt,unm1)
    nx = size(x,2);
    ny = size(y,2);
    un = zeros(nx,ny);
    for i = 1+order/2:nx-order/2
        for j = 1+order/2:ny-order/2
            if (order >= 2)
                un(i,j) = unm1(i,j) + ...
                            dt*v_t(x(i),y(j),0) + ...
                            dt^2/2*(def.c^2*(v_xx(x(i),y(j),0)+v_yy(x(i),y(j),0))+h(def,x(i),y(j),0));
%                 un(i,j) = v(x(i),y(j),dt);
            end
            if (order >= 4)
                un(i,j) = un(i,j) + ...
                            dt^3/6*(def.c^2*(v_xxt(x(i),y(j),0)+v_yyt(x(i),y(j),0))+h_t(def,x(i),y(j),0))+...
                            dt^4/24*(def.c^4*(v_4x(x(i),y(j),0)+2*v_xxyy(x(i),y(j),0)+v_4y(x(i),y(j),0))+...
                                def.c^2*(h_xx(def,x(i),y(j),0)+h_yy(def,x(i),y(j),0)));
            end
        end
    end
end

% fill in boundary conditions
% code for Dirichlet left, right, top, down BCs
% n = timestep
function un = BCs(def,x,y,order,dx,dy,un,n,dt)
    % dimensions of x and y
    mx = size(x,2);
    my = size(y,2);
    % non-ghost bounds for indices
    ia = 1+order/2;
    ja = 1+order/2;
    ib = mx-order/2;
    jb = my-order/2;
    if (order == 2)
        for i = 1:mx
            un(i,ja-1) = 2*un(i,ja)-un(i,ja+1);
            un(i,jb+1) = 2*un(i,jb)-un(i,jb-1);
            un(i,ja) = v(x(i),y(ja),n*dt);
            un(i,jb) = v(x(i),y(jb),n*dt);
        end
        for j = 1:my
            un(ia-1,j) = 2*un(ia,j)-un(ia+1,j);
            un(ib+1,j) = 2*un(ib,j)-un(ib-1,j);
            un(ia,j) = v(x(ia),y(j),n*dt);
            un(ib,j) = v(x(ib),y(j),n*dt);
        end
    elseif (order == 4)
        for j = 1:my
            un(ia,j) = v(x(ia),y(j),n*dt);
            un(ib,j) = v(x(ib),y(j),n*dt);
            % Dirichlet left BC - discrete delta fcn
            % u(1) = un(ia-1,j)
            % u(2) = un(ia-2,j)
            f = @(u) [def.c^2/dx^2*(...
                           un(ia+1,j)...
                        -2*un(ia,j)...
                          +u(1)...
                        -1/12*(...
                               un(ia+2,j)...
                            -4*un(ia+1,j)...
                            +6*un(ia,j)...
                            -4*u(1)...
                              +u(2)...
                          )...
                      ); 
                      def.c^4/dx^4*(...
                            un(ia+2,j)...
                         -4*un(ia+1,j)...
                         +6*un(ia,j)...
                         -4*u(1)...
                           +u(2)...
                      )];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0, f2-f0];
            b = -1*f0;
            u = A\b;
%             un(ia-1,j) = v(x(ia-1),y(j),n*dt);
            un(ia-1,j) = u(1);
            un(ia-2,j) = u(2);
%             un(ia-2,j) = v(x(ia-2),y(j),n*dt);
            % Dirichlet right BC - discrete delta fcn
            % u(1) = un(ib+1,j)
            % u(2) = un(ib+2,j)
            f = @(u) [def.c^2/dx^2*(...
                            un(ib-1,j)...
                         -2*un(ib,j)...
                           +u(1)...
                         -1/12*(...
                              un(ib-2,j)...
                           -4*un(ib-1,j)...
                           +6*un(ib,j)...
                           -4*u(1)...
                             +u(2)...
                         )...
                     ); 
                     def.c^4/dx^4*(...
                          un(ib-2,j)...
                       -4*un(ib-1,j)...
                       +6*un(ib,j)...
                       -4*u(1)...
                         +u(2)...
                     )];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0, f2-f0];
            b = -1*f0;
            u = A\b;
%             un(ib+1,j) = v(x(ib+1),y(j),n*dt);
%             un(ib+2,j) = v(x(ib+2),y(j),n*dt);
            un(ib+1,j) = u(1);
            un(ib+2,j) = u(2);
        end
        for i = 1:mx
            un(i,ja) = v(x(i),y(ja),n*dt);
            un(i,jb) = v(x(i),y(jb),n*dt);
            % Dirichlet top BC - discrete delta fcn
            % u(1) = un(i,jb+1)
            % u(2) = un(i,jb+2)
            f = @(u) [def.c^2/dy^2*(...
                            un(i,jb-1)...
                         -2*un(i,jb)...
                           +u(1)...
                           -1/12*(...
                                un(i,jb-2)...
                             -4*un(i,jb-1)...
                             +6*un(i,jb)...
                             -4*u(1)...
                               +u(2)...
                           )...
                      ); 
                      def.c^4/dy^4*(...
                           un(i,jb-2)...
                        -4*un(i,jb-1)...
                        +6*un(i,jb)...
                        -4*u(1)...
                          +u(2)...
                      )];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0, f2-f0];
            b = -1*f0;
            u = A\b;
%             un(i,jb+1) = v(x(i),y(jb+1),n*dt);
%             un(i,jb+2) = v(x(i),y(jb+2),n*dt);
            un(i,jb+1) = u(1);
            un(i,jb+2) = u(2);
            % Dirichlet bottom BC - discrete delta fcn
            % u(1) = un(i,ja-1)
            % u(2) = un(i,ja-2)
            f = @(u) [def.c^2/dy^2*(...
                          un(i,ja+1)...
                       -2*un(i,ja)...
                         +u(1)...
                         -1/12*(...
                              un(i,ja+2)...
                           -4*un(i,ja+1)...
                           +6*un(i,ja)...
                           -4*u(1)...
                             +u(2)...
                         )...
                      ); 
                      def.c^4/dy^4*(...
                          un(i,ja+2)...
                       -4*un(i,ja+1)...
                       +6*un(i,ja)...
                       -4*u(1)...
                         +u(2)...
                      )];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0, f2-f0];
            b = -1*f0;
            u = A\b;
%             un(i,ja-1) =  v(x(i),y(ja-1),n*dt);
%             un(i,ja-2) =  v(x(i),y(ja-2),n*dt);
            un(i,ja-1) = u(1);
            un(i,ja-2) = u(2);
        end
    end
end

% time steps for n >= 2
function unp1 = main_time_step(def,sigma_x,sigma_y,x,y,order,dt,dx,dy,unm1,un,n)
    nx = size(x,2);
    ny = size(y,2);
    unp1 = zeros(nx,ny);
    for i = 1+order/2:nx-order/2
        for j = 1+order/2:ny-order/2
            unp1(i,j) = 2*un(i,j)-unm1(i,j)+...
                sigma_x^2*(...
                       un(i-1,j)...
                    -2*un(i,j)...
                    +  un(i+1,j)...
                )+...
                sigma_y^2*(...
                       un(i,j-1)...
                    -2*un(i,j)...
                    +  un(i,j+1)...
                )+...
                dt^2*h(def,x(i),y(j),n*dt);
%             unp1(i,j) = v(x(i),y(j),(n+1)*dt);
            if (order >= 4)
                unp1(i,j) = unp1(i,j) ...
                              -sigma_x^2/12*(...
                                    un(i+2,j)...
                                 -4*un(i+1,j)...
                                 +6*un(i,j)...
                                 -4*un(i-1,j)...
                                   +un(i-2,j)...
                              ) + ...
                              -sigma_y^2/12*(...
                                    un(i,j+2)...
                                 -4*un(i,j+1)...
                                 +6*un(i,j)...
                                 -4*un(i,j-1)...
                                   +un(i,j-2)...
                              ) + ...
                              sigma_x^4/12*(...
                                    un(i+2,j)...
                                 -4*un(i+1,j)...
                                 +6*un(i,j)...
                                 -4*un(i-1,j)...
                                   +un(i-2,j)...
                              ) + ...
                              sigma_y^4/12*(...
                                    un(i,j+2)...
                                 -4*un(i,j+1)...
                                 +6*un(i,j)...
                                 -4*un(i,j-1)...
                                   +un(i,j-2)...
                              ) + ...
                              +2*def.c^4*dt^4/(12*dx^2*dy^2)*(...
                                    un(i+1,j+1)...
                                 -2*un(i+1,j)...
                                   +un(i+1,j-1)...
                                 -2*un(i,j+1)...
                                 +4*un(i,j)...
                                 -2*un(i,j-1)...
                                   +un(i-1,j+1)...
                                 -2*un(i-1,j)...
                                   +un(i-1,j-1)...
                              )...
                              +def.c^2*dt^4/12*(...
                                   h_xx(def,x(i),y(j),n*dt)...
                                  +h_yy(def,x(i),y(j),n*dt)...
                              )...
                              +dt^4/12*h_tt(def,x(i),y(j),n*dt);
                              
            end
        end
    end
end