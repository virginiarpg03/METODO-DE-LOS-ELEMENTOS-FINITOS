f = @(x,t) 0*x;

a = 0;
k = 0.1;

t0 = 0;
tf = 100;
dt = 0.00016;
Nt = (tf - t0)/dt;


c = 0;
d = 1;

Ne = 100;
h = (d-c)/Ne;
xi = c:h:d;


Calcular_Matriz_Masa_Rigidez_Lineal1D

%Euler explicito
A = M;
A0=A;
A0([1 end],:)=0;
A0(:,[1 end])=0;
A0(1,1)=1;
A0(end,end)=1;

R0=R;
R0([1 end],:)=0;
R0(:,[1 end])=0;
R0(1,1)=1;
R0(end,end)=1;

[V D]=eig(full(R0),full(M));

uhn = sin(pi*xi').^100;
plot(xi,uhn)
title(0)
pause

for n = 1:Nt
    fi = f(xi',(n-1)*dt);
    vect_b = ((1-a*dt)*M-k*dt*R)*uhn + dt*M*fi;
    vect_b([1 end])=0;

    uhn = A0\vect_b; %(A0*uh=vect_b)

    plot(xi,uhn)
    title(n*dt)
    axis([0 1 -.2 1.2])
    pause
    %Ctrol+C para parar la grafica
end
