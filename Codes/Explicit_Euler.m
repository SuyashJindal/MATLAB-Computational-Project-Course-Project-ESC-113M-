%SOLUTION USING EXPLICIT EULER METHOD
to=0; %initial time
tn=10; %final time
%hc value above which numerical stability is encountered using this method is 0.4
h=0.1;

n=(tn-to)/h;%n is the number of steps so n+1 will be the number of points

k1=2;  %reaction constant for reaction A to B
k2=0.5;  %reaction constant for reaction B to C
k3=0.3;  %reaction constant for reaction B to D

A(1)=1; %initial value of A at to=0
B(1)=0;%initial value of B at to=0
C(1)=0; %initial value of C at to=0
D(1)=0;%initial value of D at to=0
t(1)=0; %initial time      
          
for i=1:n
 A(i+1)=A(i)+h*f1(A(i),B(i),C(i),D(i));
 B(i+1)=B(i)+h*f2(A(i),B(i),C(i),D(i));
 C(i+1)=C(i)+h*f3(A(i),B(i),C(i),D(i));
 D(i+1)=D(i)+h*f4(A(i),B(i),C(i),D(i));
 t(i+1)=t(i)+h;
end

plot(t,A);
hold on;
plot(t,B);
hold on;
plot(t,C);
hold on;
plot(t,D);
hold on;
function func1= f1(A,B,C,D)  %f1= dA/dt
          k1=2;
          func1= -k1*A;
end
function func2= f2(A,B,C,D)  %f2= dB/dt
          k1=2;
          k2=0.5;
          k3=0.3;
          func2= k1*A-k2*B-k3*B;
end
function func3= f3(A,B,C,D)  %f3= dC/dt
          k2=0.5;
          func3= k2*B;
end
function func4= f4(A,B,C,D)  %f4= dD/dt
          k3=0.3;
          func4= k3*B;
end