Five iterations of Newton’s method to find the minimum of the function

syms x;

%f = 2.5*x^4-6*x^2+2*x;
%f_1st = diff(f);
%f_2nd = diff(diff(f));

fun = @(x) 2.5*x^4-6*x^2+2*x;
Dfun = @(x) 10*x^3 - 12*x + 2;
DDfun = @(x) 30*x^2 - 12;


n = 5; %number of iterations
x = zeros(n+1,1);
y = zeros(n,1);
dy = zeros(n,1);
ddy = zeros(n,1);
deltax = zeros(n,1);

x(1) = -1;

fprintf ('************************ \n');
fprintf ('NEWTON METHOD \n');
fprintf ('************************ \n');
for i =1:n
y(i) = feval(fun,x(i));
dy(i) = feval(Dfun,x(i));
ddy(i) = feval(DDfun,x(i));
deltax(i) = -dy(i)/ddy(i);
x(i+1) = x(i)+deltax(i);
fprintf('Iteration %i\n',i);
fprintf('lambda =  %i\n',x(i));
fprintf('Function = %i\n',y(i));
fprintf('change in lambda = %i\n',deltax(i));
fprintf ('************************ \n');
end

%%

n = 5; %number of iterations
x = zeros(n+1,1);
y = zeros(n,1);
dy = zeros(n,1);
ddy = zeros(n,1);
deltax = zeros(n+1,1);
x(1) = 2;
deltax(1) = 0.1;
fprintf ('************************ \n');
fprintf ('QUASI-NEWTON METHOD \n');
fprintf ('************************ \n');
for i =1:n
y(i) = feval(fun,x(i));
dy(i) = (feval(fun,x(i)+deltax(i))-feval(fun,x(i)-deltax(i)))/(2*deltax(i));
ddy(i) = (feval(fun,x(i)+deltax(i))-2*feval(fun,x(i))+feval(fun,x(i)-deltax(i)))/(deltax(i)^2);
deltax(i+1) = - dy(i)/ddy(i);
x(i+1) = x(i) + deltax(i+1);
fprintf('Iteration %i\n',i);
fprintf('lambda =  %i\n',x(i));
fprintf('Function = %i\n',y(i));
fprintf('change in lambda = %i\n',deltax(i));
fprintf ('************************ \n');
end

%%

xx = linspace(-2, 2, 50);
yy = 2.5*xx.^4-6.*xx.^2+2.*xx;
figure(1)
plot(xx, yy, '-pr')
grid
xlabel('x')
ylabel('y')
title('Function behavior')






