% analytic solution to 1-D wave eqn

function sum = waves_analytic_1d(f,g,a,b,c,x,t_f)
    L = b-a;
    A_1 = 2/L*integral(@(x) (f(x).*sin(pi*x/L)),a,b);
    B_1 = 2/(pi*c)*integral(@(x) (g(x).*sin(pi*x/L)),a,b);
    term_n = A_1.*cos(pi*c*t_f/L).*sin(pi*x/L)+B_1.*sin(pi*c*t_f/L).*sin(pi*x/L);
    sum = term_n;
    n = 2;
    while (abs(max(term_n)) > eps)
        A_n = 2/L*integral(@(x) (f(x).*sin(n*pi*x/L)),a,b);
        B_n = 2/(n*pi*c)*integral(@(x) (g(x).*sin(n*pi*x/L)),a,b);
        term_n = A_n.*cos(n*pi*c*t_f/L).*sin(n*pi*x/L)+B_n.*sin(n*pi*c*t_f/L).*sin(n*pi*x/L);
        sum = sum + term_n;
        n = n+1;
    end
end
