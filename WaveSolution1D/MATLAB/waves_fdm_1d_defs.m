% problem definitions for 1-D 

% a = left bound on x
% b = right bound on x
% c = wave number
% N = resolution in x-dimension
% t_f = final time
% u_0 = initial displacement
% v_0 = initial velocity
% left = left boundary condition: u(0,t)
% right = right boundary condition: u_x(1,t)
function [a,b,c,N,t_f,u_0,v_0,left,right] = waves_fdm_1d_defs(icase)
    switch icase
        case 1
            a = 0;
            b = 3*pi/2;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x) sin(x);
            v_0 = @(x) -sin(x);
            left = @(t) 0;
            right = @(t) 0;
        case 2
            a = 0;
            b = 3*pi/2;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x) sin(x);
            v_0 = @(x) 0;
            left = @(t) 0;
            right = @(t) 0;
        case 3
            a = 0;
            b = 1;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x) 0;
            v_0 = @(x) x;
            left = @(t) 0;
            right = @(t) 0;
        case 4
            a = 0;
            b = 1;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x) sin(x);
            v_0 = @(x) -sin(x);
            left = @(t) t;
            right = @(t) 0;
        case 5
            a = 0;
            b = 1;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x) sin(x);
            v_0 = @(x) -sin(x);
            left = @(t) 0;
            right = @(t) 0;
    end
end
