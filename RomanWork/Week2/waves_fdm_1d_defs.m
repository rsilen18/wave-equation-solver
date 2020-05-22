% problem definitions for 1-D 

% a = left bound on x
% b = right bound on x
% u_0 = initial displacement
% v_0 = initial velocity
% left = left boundary condition
function [a,b,u_0,v_0,left] = waves_fdm_1d_defs(icase)
    switch icase
        case 1
            a = 0;
            b = 3*pi/2;
            u_0 = @(x) sin(x);
            v_0 = @(x) -sin(x);
            left = @(t) 0;
        case 2
            a = 0;
            b = 3*pi/2;
            u_0 = @(x) sin(x);
            v_0 = @(x) 0;
            left = @(t) 0;
        case 3
            a = 0;
            b = 1;
            u_0 = @(x) 0;
            v_0 = @(x) x;
            left = @(t) 0;
        case 4
            a = 0;
            b = 1;
            u_0 = @(x) sin(x);
            v_0 = @(x) -sin(x);
            left = @(t) t;
    end
end