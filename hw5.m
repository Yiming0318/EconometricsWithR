%HW4
clear all;close all;clc;
%problem1
x = @(t) 1/2*exp(-t/2)+2*t*exp(-t/2);
[t_max x_tmax] = fminbnd(@(t) -x(t),0,5);
A1 = [t_max abs(x_tmax)];
save('A1.dat','A1','-ascii');

%problem2
a=0; b=5;
c=(-1+sqrt(5))/2;
t1 =c*a + (1-c)*b;
t2 = (1-c)*a+c*b;
f1 = -(1/2*exp(-t1/2)+2*t1*exp(-t1/2));
f2 = -(1/2*exp(-t2/2)+2*t2*exp(-t2/2));
for j = 1:100
    if f1 < f2
        b = t2; t2=t1;f2=f1;
        t1 = c*a+(1-c)*b;
        f1 = -(1/2*exp(-t1/2)+2*t1*exp(-t1/2));
    else
        a = t1;t1=t2;f1=f2;
        t2 = (1-c)*a +c*b;
        f2 = -(1/2*exp(-t2/2)+2*t2*exp(-t2/2));
    end
    if (b - a)<10^(-3);
        break
    end
end
A2=[a,b];
save('A2.dat','A2','-ascii');

%problem3
f = @(x)2*exp(-x/2)+(2*x+1/2)*(-1/2)*exp(-1/2*x);
fprime = @(x)-exp(-1/2*x)-(exp(-1/2*x)+(x+1/4)*(-1/2)*exp(-1/2*x));
xk = 0;
tol = 10^(-3);
maxIter = 10000;
[x,iter] = Newton(f,fprime,xk,tol,maxIter);
save('A3.dat','x','-ascii');
save('A4.dat','iter','-ascii');


%prbolem4
syms x y real
A = [x; y];
f = @(A)(1-A(1))^2+2*(A(2)-A(1)^2)^2;
A5 = fminsearch(f, [1;10]);
save('A5.dat','A5','-ascii');

%problem 5
fgradxy = @(x,y) [-2 *(1-x)- 8 * x *(y - x^2); 4 * (y - x^2)];
fgrad = @(p) fgradxy(p(1),p(2));
p = [1; 10];
step = 0;
while step < 10000 
    step = step + 1;
    grad = fgrad(p); 
    phi = @ (t) p - t.*grad; 
    f_of_phi = @ (t) f(phi(t)); 
    tmin = fminbnd(f_of_phi,0,1); 
    p = phi(tmin);
    if norm(grad,inf) < 10^(-4)
       break 
    end
end
A7 = step-1;
save('A6.dat','p','-ascii');
save('A7.dat','A7','-ascii');

%problem6
tstep = 0.01;
while step < 10000 
    step = step + 1;
    grad = fgrad(p); 
    p = p-tstep*grad;
    if norm(grad,inf) < 10^(-4)
       break 
    end
end
fminsearch(f,[1;10]);
A8 = step -1;
save('A8.dat','A8','-ascii');

tstep = 0.02;
while step < 10000 
    step = step + 1;
    grad = fgrad(p); 
    p = p-tstep*grad;
    if norm(grad,inf) < 10^(-4)
       break 
    end
end
A9 = step - 1;
save('A9.dat','A9','-ascii');

tstep = 0.025;
while step < 10000 
    step = step + 1;
    grad = fgrad(p); 
    p = p-tstep*grad;
    if norm(grad,inf) < 10^(-4)
       break 
    end
end
A10 = step - 1;
save('A10.dat','A10','-ascii');



