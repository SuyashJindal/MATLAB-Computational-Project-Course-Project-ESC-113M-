k1 = 2;
k2 = 0.5;
k3 = 0.3;

y0 = [1;0;0;0];
ivp_implicit_eulersys(0,y0,0.1,100)
solution

function yout = ivp_implicit_eulersys(x0,y0,h,n_out)
xout = zeros(n_out+1,1);
xout(1) =x0;
yout = zeros(length(y0),n_out+1);
%yout = zeros(n_out+1,length(y0));
%yout(1,:) = y0;
yout(:,1) = y0;
x = x0;
y = y0;
tol = 0.01;
for j = 2:n_out+1
    yold = y;
    R = getR(y,yold,h);
    count = 1;
    while norm(R)>tol
        J = getJ(y,h);
        del = -J\R;
        y = y+del;
        R = getR(y,yold,h);
        count = count+1;
        if count>1000
            fprintf('Failed')
            return
        end
    end
    x= x+h;
    xout(j) = x;
    %yout(j,:) = y;
    yout(:,j) = y;
end
plot(xout,yout)
title('Implicit Method');
end
function R = getR(y,yold,h)
k1 = 2;
k2 = 0.5;
k3 = 0.3;

R  = zeros(4,1);
R(1)= y(1)-yold(1)+h*(k1*y(1));
R(2) = y(2)-yold(2)-h*(k1*y(1)-k2*y(2)-k3*y(2));
R(3)= y(3)-yold(3)-h*(k2*y(2));
R(4) = y(4)-yold(4)-h*(k3*y(2));
end
function J = getJ(y,h)
J = zeros(4);
k1 = 2;
k2 = 0.5;
k3 = 0.3;
J(1,1) = 1+h*(k1);
J(1,2) = 0;
J(1,3) = 0;
J(1,4) = 0;
J(2,1) = h*k1;
J(2,2) = 1 -h*(-k2-k3);
J(2,3) = 0;
J(2,4) = 0;
J(3,1) = 0;

J(3,2) =  -h*(k2);
J(3,3) = 1;
J(3,4) = 0;
J(4,1) = 0;
J(4,2) = -h*(k3);
J(4,3) = 0;
J(4,4) = 1;
end
function solution
[t,y] = ode45(@getf,[0,10],[1,0,0,0])

plot(t,y(:,1),'r')
figure(2)
plot(t,y(:,2),'g')
figure(3)
plot(t,y(:,3),'b')
figure(4)
plot(t,y(:,4),'o')
figure(5)
plot(t,y,'LineWidth',5)
xlabel('Time')
ylabel('Concentrations')



end
function f = getf(t,y)
k1 = 2;
k2 = 0.5;
k3 = 0.3;
f(1,1) = -k1*y(1);
f(2,1) = k1*y(1)-k2*y(2)-k3*y(2);
f(3,1) = k2*y(2);
f(4,1) = k3*y(2);
end






