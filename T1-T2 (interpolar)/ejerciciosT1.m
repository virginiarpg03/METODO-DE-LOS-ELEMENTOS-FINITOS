%% EJERCICIOS HOJA 1:


%% ejercicio1

f = @(x) (3 + x).*(cos(pi.*x./4)).^2;

xi=[0 1 3];
fi= f(xi);

x=[2, 2.4 3.5 4];

p = interpol_lagrange(xi,fi,x) %valores en esos puntos

% si quiero dibujarlo :

x=0:0.001:4;
p = interpol_lagrange(xi,fi,x);
hold on
plot(x,f(x),'r') %funcion real
plot(x,p, 'b')
hold off

%% ejercicio 2  hacerlo en casa con la formula que no me acuerdo

%error lagrange= max(derivada n)*w/facorial(n)



%% ejercicio 3 
clear 
% a) interpolación lineal: con 2 nodos, (1o y 3o) 
xi= [9 10];
fi= [2.1972  2.3026];

x= [9.2]; %punto aproximado

p = interpol_lagrange(xi,fi,x)

% dibujamos para ver que tal la aproximacion
x = 8.5:0.001:10.5;
p = interpol_lagrange(xi,fi,x);
hold on
plot(x,log(x),'r')
plot(x,p,'b')
hold off


%cota de error a mano con la fórmula

% b) interpolación cuadrática: con los 3 nodos

xi= [9 9.5 10];
fi= [2.1972 2.2513 2.3026];

x= [9.2]; %punto aproximado

p = interpol_lagrange(xi,fi,x)

% dibujamos para ver que tal la aproximacion
x = 8.5:0.001:10.5;
p = interpol_lagrange(xi,fi,x);
figure
hold on
plot(x,log(x),'r')
plot(x,p,'g')
hold off

%cota de error a mano con la fórmula


%% ejercicio 4
f = @(x) 2.^x;
xi = [-1 0 1];
fi =f(xi);

x=[-1.7];

p=interpol_lagrange(xi,fi,x)

real= f(x)      %valor real: 

%cota del error

% dibujamos para ver que tal la aproximacion
x = -1:0.001:1;
p = interpol_lagrange(xi,fi,x);
figure
hold on
plot(x,log(x),'r')
plot(x,p,'b')
hold off


%% ejercicio 5

f=@(x) x.^5 + 3*x.^12;
xi=0:1:20;
fi=f(xi);
%revisa
p = interpol_lagrange(xi,fi,xi);

x=1:0.001:20;
%p = interpol_lagrange(xi,fi,x);
hold on
plot(x,f(x),'r')
plot(xi,p,'b')
hold off

%% ejercicio 6

f=@(x) 1./(1+x.^2);
xi= linspace(-5,5,15);
fi=f(xi);

%revisa
p = interpol_lagrange(xi,fi,xi);

x= -5:0.001:5;
p = interpol_lagrange(xi,fi,x);
hold on
plot(x,f(x),'r')
plot(x,p,'b')
hold off

%% ejercicio 7

f=@(x) 1./(1+x.^2);
n=14;

i=0:1:n;

xj=cos((2*i + 1)*pi./(2*n +2));

%cambio de variable:
xi=xj*5;
fi=f(xi);
x= -5:0.001:5;

p=interpol_lagrange(xi,fi,x);

figure
hold on
plot(x,f(x),'r')
plot(x,p,'b')
plot(xi,fi,'o')
hold off



%% ejercicio 8

% a) f(0.8)-f(0.8)
%lineal: solo dos puntos

f=@(x) sin(pi*x);

xi=[0.5 1];

fi=f(xi);

p=interpol_lagrange(xi,fi,0.8);

error_lineal = abs(p-f(0.8))

%cuadrático 3 puntos

xi = [0 0.5 1];
fi = f(xi);

p=interpol_lagrange(xi,fi,0.8);
error_cuadratico = abs(p-f(0.8))

%% ejercicio 9

f=@(x) 1./(1+x.^2);

% grado 1: en cada elemento 2 nodos
figure ('Name','aproximacion lineal')
Ne=7;                                 %numero elementos
Nn= Ne+1;                             %numero nodos
malla= linspace(-5,5,Nn);
for i=1:7
    xi=[malla(i) malla(i+1)];
    fi=f(xi);
    x=malla(i):0.001:malla(i+1);
    p=interpol_lagrange(xi,fi,x);

    hold on
    plot (x,f(x),'k') %real
    plot (x,p)
    hold off
end

% grado 2: en cada elemento 3 nodos
figure ('Name','aproximación cuadrática')
m=2;
a=-5;
h=10/7;         %distancia entre elementos
for i=1:7
    b=a+h;
    xi=linspace(a,b,m+1);
    fi=f(xi);
    x=a:0.001:b;
    p=interpol_lagrange(xi,fi,x);

    hold on
    plot (x,f(x),'k') %real
    plot (x,p)
    hold off
    a=b %avanzar
end


%grado 3: en cada elemento 4 nodos

figure ('Name','aproximación cúbica')
m=3;
a=-5;
for i=1:7
    b=a+h;
    xi=linspace(a,b,m+1);
    fi=f(xi);
    x=a:0.001:b;
    p=interpol_lagrange(xi,fi,x);

    hold on
    plot (x,f(x),'k') %real
    plot (x,p)
    hold off
    a=b %avanzar
end


%% ejercicio 10
clear all
close all

f=@(x) 1./(1+x.^2);
a=-5;
b=5;

%a) h=5/6 m=1
m=1
Ne=12;
elementos_finitos_1D(f,a,b,Ne,m);
title('a)')

%c) h=10/3 m=4
m=4
Ne=3;
elementos_finitos_1D(f,a,b,Ne,m);
title('c)')

% random
m=14
Ne=10;
elementos_finitos_1D(f,a,b,Ne,m);
title('uwu')


%% ejercicio 11




















