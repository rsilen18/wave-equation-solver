% code to demonstrate CFL condition

[xi,sigma] = meshgrid(-pi:pi/100:pi,-1.1:0.01:1.1);
% discr = (sigma.*cos(xi)).^2-(sigma.^2-2).*cos(xi)+(sigma.^2-2);
% A_plus = sigma.^2.*cos(xi)+sigma.*sqrt(discr);
% A_minus = sigma.^2.*cos(xi)-sigma.*sqrt(discr);
discr = (2*sigma.^2-2*sigma.^2.*cos(xi)-2).^2-4;
A_plus = (-2*sigma.^2+2*sigma.^2.*cos(xi)+2+sqrt(discr))/2;
A_minus = (-2*sigma.^2+2*sigma.^2.*cos(xi)+2-sqrt(discr))/2;


% xi = linspace(-pi,pi,101);
% sigma = linspace(-1.1,1.1,101);
% A_plus = zeros(101);
% A_minus = zeros(101);
% for i = 1:101
%     for j = 1:101
%         discr = (2*sigma(j)^2-2*sigma(j)^2*cos(xi(i))-2)^2-4;
%         A_plus(i,j) = (-2*sigma(j)^2+2*sigma(j)^2*cos(xi(i))+2+sqrt(discr))/2;
%         A_minus(i,j) = (-2*sigma(j)^2+2*sigma(j)^2*cos(xi(i))+2-sqrt(discr))/2;
%     end
% end

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

% figure(3)
% A = abs(A_plus.*A_minus);
% surf(xi,sigma,A);
% shading interp
% xlabel("xi");
% ylabel("sigma");
% zlabel("|A|");