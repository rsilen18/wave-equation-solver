% problem definitions for 2-D 

% a_x = lower bound on x
% b_x = upper bound on x
% a_y = lower bound on y
% b_y = upper bound on y
% c = wave number
% N = resolution in x-dimension
% t_f = final time
% u_0 = initial displacement: u(x,y,0)
% v_0 = initial velocity: u_t(x,y,0)
% left = left boundary condition: u(a_x,y,t)
% right = right boundary condition: u(b_x,y,t)
% bottom = bottom BC: u(x,a_y,t)
% top = top BC: u(x,b_y,t)
function [a_x,b_x,a_y,b_y,c,N,t_f,u_0,v_0,left,right,bottom,top] = waves_fdm_2d_defs(icase)
    switch icase
        case 1
            a_x = 0;
            b_x = 1;
            a_y = 0;
            b_y = 1;
            c = 1;
            N = 100;
            t_f = 10;
            u_0 = @(x,y) x*y*(1-x)*(1-y);
            v_0 = @(x,y) sin(2*pi*x);
            left = @(y,t) 0;
            right = @(y,t) 0;
            bottom = @(x,t) 0;
            top = @(x,t) 0;
    end
end