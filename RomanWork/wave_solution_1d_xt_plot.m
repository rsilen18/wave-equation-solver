% set-up
N = input("Enter N: ");     % resolution in x-direction
c = input("Enter c: ");     % wave speed
t_f = input("Enter final time: ");
a = 0;  % left bound on x
b = 3*pi/2;  % right bound on x
x = linspace(a,b,N+1);
dx = (b-a)/N;
dt = .01;
sigma = c*dt/dx
assert(x(2)-x(1)==dx);

% evaluation
u = zeros(N+2,1);
% step 1 - initial condition
for j = 1:N+1
   u(j) = u_0(x(j)); 
end
% step 2 - initial ghost point
u(N+2) = u(N);
% step 3 - first time step
new = zeros(N+2,1);
new(1) = 0;
for j = 2:N+1
   new(j) = (1-sigma^2)*u(j)+dt*v_0(x(j))+sigma^2/2*(u(j-1)+u(j+1)); 
end
% step 4 - first time step ghost point
new(N+2) = new(N);
u = [u,new];
% step 5 - n = time, j = space
for n = 2:(t_f/dt+1)
   new = zeros(N+2,1);
   % step 6
   for j = 2:N+1
        new(j) = 2*u(j,n)-u(j,n-1)+sigma^2*(u(j+1,n)-2*u(j,n)+u(j-1,n));
   end
   new(N+2) = new(N);
   u = [u,new];
end

% plotting

% 3-d plot
figure(1);
surf(u);
shading interp
xlabel("x");
ylabel("t");
zlabel("u");

% 2-d animation
for t = 1:(t_f/dt+1)
   figure(2);
   plot(x,u(1:N+1,t)); 
   pause(.01);
end


% given functions
function value = u_0(x)
    value = sin(x);
end

function value =  v_0(x)
    value = -sin(x);
end
