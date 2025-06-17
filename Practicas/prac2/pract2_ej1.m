%Ejercicio 1
f=@(x) 2+sin(pi*x);
a=1;   %coef de u
k=0.005;   %coef de u segunda

c=0;   %limite inferior intervalo
d=2;   %limite superior intervalo

Ne=100;
h=(d-c)/Ne;
xi=c:h:d;

calcular_matriz_masa_rigidez_lineal
A=a*M+k*R;
A0=A;
A0([1 end],:)=0;
A0(:,[1 end])=0;
A0(1,1)=1;
A0(end,end)=1;

fi=f(xi');

gi=0*xi';
vect_b= M*fi-A*gi;
vect_b([1 end])=0;


wh=A0\vect_b;     
uh=wh+gi;
plot(xi,uh,'o-')    %dibuja la solución aproximada

