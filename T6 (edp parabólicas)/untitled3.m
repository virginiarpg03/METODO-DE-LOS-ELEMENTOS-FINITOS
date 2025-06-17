k=1e-2; a= 0; f= @(x,t) 50*exp(-x.^2/0.01); %Parámetros EDP
c = -1; d = 4; %Intervalo espacial
t0 = 0; tf = 100; %Intervalo temporal 

%Partición temporal para resolución por PVI
dt = 0.000016;
Nt = (tf-t0)/dt;

%Partición espacial para resolución por MEF
Ne = 1000;
h = (d-c)/Ne;
xi = c:h:d;

%RESOLUCIÓN
Calcular_Matriz_Masa_Rigidez_Lineal1D
%Valor de la solución para t0
uhn = 0*xi'+20; %Para todos los puntos que no pertenecen a la frontera, whn=uhn
plot(xi,uhn)
title(0)
pause(0.05)

%Hago el cálculo de A y A0 aquí fuera porque no varía
A =(1+0.5*a*dt)*M + 0.5*k*dt*R; A0 = A; 
A0([1 end],:) = 0;
A0(:,[1 end]) = 0;
A0(1,1) = 1;
A0(end,end) = 1;

%Defino mi función auxiliar
gi=0*xi';
gi(1)=5;
gi(end)=30;
whn=uhn+gi;

figure
fiold=f(xi',0);

%Para cada instante de tiempo tengo un problema elíptico
for n=1:Nt
    fi = f(xi',n*dt);

    vect_b = ((1-a*dt*0.5)*M-(k*dt*0.5)*R)*whn + dt*0.5*M(fi+fiold)-A*gi;
    vect_b0=vect_b;
    vect_b0([1 end])=0; %Condición de Dirichlet nula

    whn = A0\vect_b0; 
    uhn= whn + gi;
     
    plot(xi,uhn)
    title(n*dt)
    axis([-1 4 -5 35])
    pause(0.05)

    fiold=fi;
end