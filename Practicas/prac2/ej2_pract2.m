f=@(x) 2+sin(pi*x);
a=1;   %coef de u
k=0.005;   %coef de u segunda

c=0;   %limite inferior intervalo
d=2;   %limite superior intervalo

Ne=100;
h=(d-c)/Ne;
xi=c:h:d;

Calcular_Matriz_Masa_Rigidez_Lineal1D
A=a*M+k*R;
A0=A;

C=zeros(Ne+1,Ne+1);
C(end,end)=1;

fi=f(xi');

gi=0*xi';
vect_b= M*fi-A*gi;

uh=(A0+C)\vect_b; 
plot(xi,uh,'-o')  

uh(end)  %solución evaluada en 2 (uh(2))