function u = genSolverA(inter,N,tf,sigma,c,f,g,plot_)
    x = linspace(inter(1), inter(2), N+1);
    dx = x(2)-x(1);
    M = ceil(tf/(sigma*dx/c));
    dt = tf/M;
    sigma = c*dt/(x(2)-x(1));
    
    [u1, u2] = genICA(x,f,g,sigma,dt);
    if (plot_)
        plotter(x,u1,tf/M);
        plotter(x,u2,tf/M);
    end
    
    
    unm1 = u1; 
    un = u2;
    
    for n = 3:M
        unp1 = genStep(un,unm1);
        if (plot_)
            plotter(x,unp1,tf/M);
        end
        unm1 = un;
        un = unp1;
    end
    u = un;
end

function unp1 = genStep(un,unm1,sigma)
    unp1 = un;
    un(1) = 0;
    for j = 2:N
        unp1(j) = 2*un(j)-unm1(j)+sigma^2*(un(j+1)-2*un(j)+un(j-1));
    end
    unp1(N+1) = (4*unp1(N)-unp1(N-1))/3;
end

function plotter(x,u,wait)
    figure(1)
    plot(x,u)
    xlabel('x (Spacial Dimension)')
    ylabel('u (Functional Value)')
    pause(wait);
end