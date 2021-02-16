% plot symbols and errors for 1st, 2nd, and 4th order difference methods

delta_x = 1;
xi = linspace(-pi,pi,100);
k = xi/delta_x;

symbol_2 = 2*cos(xi)-2;
symbol_4 = (2*cos(xi)-2)-(1/12)*(2*cos(xi)-2).^2;
d_xx = -k.^2;

figure(1)
plot(xi,d_xx);
hold on
plot(xi,symbol_2);
plot(xi,symbol_4);
legend("d_{xx}","Re(S_2)","Re(S_4)");
title("Re(S) vs xi plot");
xlabel("xi");
ylabel("Re(S)");
hold off;

figure(2)
plot(xi,abs(d_xx-symbol_2));
hold on
plot(xi,abs(d_xx-symbol_4));
legend("Re(S_2)-(-k^2)","Re(S_4)-(-k^2)");
xlabel("xi");
ylabel("Error");
title("Error plot: d_xx");
hold off

d_x = imag(1i*k);
symbol_1 = imag(cos(xi)+1i*sin(xi)-1);

figure(3)
plot(xi,d_x);
hold on
plot(xi,symbol_1);
legend("d_x","Imag(S_1)");
title("Imag(S) vs xi plot");
xlabel("xi");
ylabel("Re(S)");
hold off;

figure(4)
plot(xi,abs(d_x-symbol_1));
xlabel("xi");
ylabel("Error");
title("Error plot: d_x");
