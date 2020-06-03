% code to demonstrate CFL condition

% if false, look at 4th order stability
order2 = false;

[xi,sigma] = meshgrid(-pi:pi/100:pi,-1.1:0.01:1.1);
% term for D+D-D+D-u
o4_term = 1/24*sigma.^2.*(-1*exp(2i*xi)+16*exp(1i*xi)-30+16*exp(-1i*xi)-exp(-2i*xi));
% correction term
corr = 1/24*sigma.^4.*(exp(2i*xi)-4*exp(1i*xi)+6-4*exp(-1i*xi)+exp(-2i*xi));

if (order2 == true)
    b = 1+sigma.^2.*(cos(xi)-1);
else
    b = 1+o4_term+corr;
end
A_plus = b+sqrt(b.^2-1);
A_minus = b-sqrt(b.^2-1);


% get magnitudes
A_plus = sqrt(real(A_plus).^2+imag(A_plus).^2);
A_minus = sqrt(real(A_minus).^2+imag(A_minus).^2);

figure(1)
surf(xi,sigma,A_plus);
shading interp
colorbar
xlabel("xi");
ylabel("sigma");
zlabel("|A_+|");

figure(2)
surf(xi,sigma,A_minus);
shading interp
colorbar
xlabel("xi");
ylabel("sigma");
zlabel("|A_-|");
