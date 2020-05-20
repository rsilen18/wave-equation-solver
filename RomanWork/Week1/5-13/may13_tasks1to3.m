N = input("Enter N: ");     % resolution in x-direction
k = input("Enter k: ");     % wave number
c = input("Enter c: ");     % wave speed
t_f = input("Enter final time: ");
a = 0;  % left bound on x
b = 1;  % right bound on x
x = linspace(a,b,N+1);
dx = (b-a)/N;
assert(x(2)-x(1)==dx);
dt = t_f/100;
for n = 1:100
   t = n*dt; 
   u = cos(k*(x-c*t));
   plot(x,u);
   xlim([a b]);
   ylim([-1 1]);
   str = sprintf('Left to right wave (t = %.1f)', t);
   title(str);
   xlabel("x");
   ylabel("u(x,t)");
   pause
end

for n = 1:100
   t = n*dt; 
   u = cos(k*(x+c*t));
   plot(x,u);
   xlim([a b]);
   ylim([-1 1]);
   str = sprintf('Right to left wave (t = %.1f)', t);
   title(str);
   xlabel("x");
   ylabel("u(x,t)");
   pause
end

for n = 1:100
   t = n*dt; 
   u = cos(k*(x-c*t))+cos(k*(x+c*t));
   plot(x,u);
   xlim([a b]);
   ylim([-2 2]);
   str = sprintf('Standing wave (t = %.1f)', t);
   title(str);
   xlabel("x");
   ylabel("u(x,t)");
   pause
end