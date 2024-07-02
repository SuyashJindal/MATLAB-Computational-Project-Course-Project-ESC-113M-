%SOLUTION USING RK4 METHOD

T1=0;  %T1=initial time
Tend=10;   %Tend=final time

%{
only those values of h are allowed for which N i.e. number of data points 
come out to be integers,
such values are 0.1,0.2,0.4,0.5,1,2,etc.
hc value is the value after which numerical unstability occurs for allowed
values of h,
hc value for the question using RK4 method is 1. 
hc=1.
%}
h=0.4;   

N= 1+ (10/h);  %N= Number of data points
time= linspace(T1,Tend,N); %time vector for particular question having N data points
k1=2;  %reaction constant for reaction A to B
k2=0.5;  %reaction constant for reaction B to C
k3=0.3;  %reaction constant for reaction B to D
Ao=1;   %initial value of A at T1=0
Bo=0;  %initial value of B at T1=0
Co=0;  %initial value of C at T1=0
Do=0;  %initial value of D at T1=0
Yold= [Ao; Bo; Co; Do];  %vector created to denote values of A,B,C,D at some ith time (1<=i<=N)
Ynew= zeros(4,1); %vector created to denote values of A,B,C,D at (i+1)th time
A= zeros(1,N); %vector created to store values of A 
B= zeros(1,N); %vector created to store values of B
C= zeros(1,N); %vector created to store values of C
D= zeros(1,N); %vector created to store values of D
A(1)=Ao;
B(1)=Bo;
C(1)=Co;
D(1)=Do;

for i= 2:N
    T= (i-1)*h;  %T=time for which concentration values of A,B,C,D is being calculated
    a= Yold(1);
    b= Yold(2);
    c= Yold(3);
    d= Yold(4);

    a1= h*f1(a,b,c,d,T);
    b1= h*f2(a,b,c,d,T);
    c1= h*f3(a,b,c,d,T);
    d1= h*f4(a,b,c,d,T);

    a2= h*f1(a+0.5*a1,b+0.5*b1,c+0.5*c1,d+0.5*d1,T+0.5*h);
    b2= h*f2(a+0.5*a1,b+0.5*b1,c+0.5*c1,d+0.5*d1,T+0.5*h);
    c2= h*f3(a+0.5*a1,b+0.5*b1,c+0.5*c1,d+0.5*d1,T+0.5*h);
    d2= h*f4(a+0.5*a1,b+0.5*b1,c+0.5*c1,d+0.5*d1,T+0.5*h);
    
    a3= h*f1(a+0.5*a2,b+0.5*b2,c+0.5*c2,d+0.5*d2,T+0.5*h);
    b3= h*f2(a+0.5*a2,b+0.5*b2,c+0.5*c2,d+0.5*d2,T+0.5*h);
    c3= h*f3(a+0.5*a2,b+0.5*b2,c+0.5*c2,d+0.5*d2,T+0.5*h);
    d3= h*f4(a+0.5*a2,b+0.5*b2,c+0.5*c2,d+0.5*d2,T+0.5*h);

    a4= h*f1(a+a3,b+b3,c+c3,d+d3,T+h);
    b4= h*f2(a+a3,b+b3,c+c3,d+d3,T+h);
    c4= h*f3(a+a3,b+b3,c+c3,d+d3,T+h);
    d4= h*f4(a+a3,b+b3,c+c3,d+d3,T+h);
    
    anew= a+ (h/6)*(a1+2*a2+2*a3+a4);
    bnew= b+ (h/6)*(b1+2*b2+2*b3+b4);
    cnew= c+ (h/6)*(c1+2*c2+2*c3+c4);
    dnew= d+ (h/6)*(d1+2*d2+2*d3+d4);

    Ynew(1)=anew;
    Ynew(2)=bnew;
    Ynew(3)=cnew;
    Ynew(4)=dnew;

    A(1,i)= Ynew(1);  
    B(1,i)= Ynew(2);
    C(1,i)= Ynew(3);
    D(1,i)= Ynew(4);

    Yold= Ynew;  %updating value of Yold for next iteration (since i will change to i+1)
  
end

figure(1)
plot(time,A)
legend('A')
hold on
plot(time,B,'DisplayName','B')
hold on
plot(time,C,'DisplayName','C')
hold on
plot(time,D,'DisplayName','D')
hold off

function func1= f1(A,B,C,D,T)  %f1= dA/dt
          k1=2;
          func1= -k1*A;
end
function func2= f2(A,B,C,D,T)  %f2= dB/dt
          k1=2;
          k2=0.5;
          k3=0.3;
          func2= k1*A-k2*B-k3*B;
end
function func3= f3(A,B,C,D,T)  %f3= dC/dt
          k2=0.5;
          func3= k2*B;
end
function func4= f4(A,B,C,D,T)  %f4= dD/dt
          k3=0.3;
          func4= k3*B;
end
