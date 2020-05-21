function [u1,u2] = genICA(x,f,g,sigma,dt)
    u1 = f(x);
    u2 = u1;
    for j = 2:(size(x,2)-1)
        u2(j) = (1-sigma^2)*u1(j) + dt*g(x(j)) + sigma^2/2*(u1(j-1)+u1(j));
    end
    u2(N+1) = (4*u2(N)-u2(N-1))/3;
end