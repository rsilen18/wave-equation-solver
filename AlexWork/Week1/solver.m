function u = solver(u_curr,u_old,t,dt,x,dx,c)
    if (t == 0)
        u = t_0(x);
    elseif (t == dt)
        sigma = c*dt/dx;
        u = t_1(u_curr, sigma, dt, dx);
    else
        sigma = c*dt/dx;
        u = t_n(u_curr,u_old,sigma);
    end
end

function u = t_0(x)
    u = f(x);
    u(size(x,2)) = u(size(x,2)-2);
end

function u = t_1(u_curr, sigma, dt, dx)
  u = u_curr;
  u(1) = 0;
  for i = 2:(size(u,2)-1)
      u(i) = (1-sigma^2)*u_curr(i) + dt*(v(i*dx)) +sigma^2/2*(u_curr(i-1) + u_curr(i+1));
  end
  u(size(u,2)) = u(size(u,2)-2);
end

function u = t_n(u_curr,u_old,sigma)
    u = u_curr;
    u(1) = 0;
    for i = 2:(size(u,2)-1)
      u(i) = 2*u_curr(i) - u_old(i) + sigma^2*(u_curr(i-1) - 2*u_curr(i) + u_curr(i+1));
    end
    u(size(u,2)) = u(size(u,2)-2);
end

function y = f(x)
    y = sin(3*pi/2*x);
end

function y = v(x)
    y = -sin(3*pi/2*x);
end