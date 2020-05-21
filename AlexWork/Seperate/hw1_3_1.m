
a = 0;
b = 0;
M = 10;
nu = 1/2;
t_del = 0.02;
t_max = 0.1;

x = ones(M+1,1);
for i = 1:size(x)
    x(i) = 1/M*(i-1);
end
x
y = make(x);
figure(1)
plot(x,y)
hold on

t = t_del;
y_new = y;
y_old = y;

for i = 1:size(y,1)
    if (i == 1)
        y(i) = 0;
    elseif (i == size(y_new,1))
        y(i) = 0;
    else
        y(i) = y_old(i) + nu*t_del*M*M*(y_old(i+1)-2*y_old(i)+y_old(i-1));
    end
end


i = 0;
while (t <= t_max)
    y
    plot(x,y)
    for i = 1:size(y_new,1)
        if (i == 1)
            y_new(i) = 0;
        elseif (i == size(y_new,1))
            y_new(i) = 0;
        else
            y_new(i) = y_old(i) + 2*nu*t_del*M*M*(y(i+1)-2*y(i)+y(i-1));
        end
    end
    y_old = y;
    y = y_new;
    t = t + t_del;
end


function out = make(x)
    out = ones(size(x,1),1);
    for i = 1:size(x)
        out(i) = f(x(i));
    end
end

function y = f(x)
    y = sin(2*pi*x);
end