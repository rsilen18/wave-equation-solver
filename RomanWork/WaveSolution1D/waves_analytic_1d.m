% analytic solution to 1-D wave eqn

% d'Alembert's formula
% u = u(x,t_f)
function u = waves_analytic_1d(f,g,c,x,t_f)
    if (size(x,2) > 1)
        u = zeros(1,size(x,2));
        u(1) = f(0);
        for i = 2:size(x,2)
            u(i) = 0.5*(f(x(i)-c*t_f)+f(x(i)+c*t_f))+1/(2*c)*integral(g,x(i)-c*t_f,x(i)+c*t_f);
        end
    else
        u = 0.5*(f(x-c*t_f)+f(x+c*t_f))+1/(2*c)*integral(g,x-c*t_f,x+c*t_f);
    end
end
