% Visualization of analytical solutions to 1-D Wave Equation using MatLab
% Author: Alexandre Ait-Ettajer
% Date: May 19, 2020 

% User Input of Parameters for wave equation solutions
N = input('Resolution of x-space (0,1) (N): ');
M = input('Total Time Steps (M): ');
k = input('Wave Parameter (k): ');
c = input('Wave Speed (c): ');
t_f = input('Final Time (t_f): ');
a = 0;
b = 1;

% Define x and t interval
x = linspace(a,b,(N));
dx = (b-a)/(N-1);
dt = t_f/M;
assert(x(2)-x(1) == dx);

% Loop over 
for n = 1:M
   t = n*dt;
   u1 = f1(x,t,c,k);
   u2 = f2(x,t,c,k);
   u3 = f3(x,t,c,k);

   % Figure of Solution 1 (Forward Propagating Wave)   
   figure(1)
   plot(x,u1);
   xlim([a b])
   ylim([-1 1])
   title('Wave Equation Solution: Backward Propagation')
   xlabel('x')
   ylabel('u')
   
   % Figure of Solution 2 (Backward Propagating Wave)
   figure(1)
   plot(x,u2);
   xlim([a b])
   ylim([-1 1])
   title('Wave Equation Solution: Forward Propagation')
   xlabel('x')
   ylabel('u')   
   
   % Figure of Solution 3 (Standing Wave)  
   figure(1)
   plot(x,u3);
   xlim([a b])
   ylim([-1 1])
   title('Wave Equation Solution: Standing Wave')
   xlabel('x')
   ylabel('u')

   pause(t_f/M);
end

% Solution 1 to the 1-D wave eqaution
function u = f1(x,t,c,k)
    u = 0.5*cos(2*pi*k*(x+c*t));
end

% Solution 2 to the 1-D wave eqaution
function u = f2(x,t,c,k)
    u = 0.5*cos(2*pi*k*(x-c*t));
end

% Solution 3 to the 1-D wave eqaution (Superposition of Sol. 1 and 2)
function u = f3(x,t,c,k)
    u = f1(x,t,c,k) + f2(x,t,c,k);
end