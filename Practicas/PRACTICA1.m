%% PRÁCTICA 1 


%% ejercicio 1
%%1
% función que, dado el numero de elementos Ne, el grado m, el intervalo (a,b) y la
% función f, dibuja la interpolación y la función real

function [p] = elementos_finitos_1D (f,a,b,Ne,m)

h = (b-a)/Ne; %distancia entre elementos


for i=1:Ne
    c=a+h;
    xi=linspace(a,c,m+1);
    fi=f(xi);
    x=a:0.001:c;
    p=interpol_lagrange(xi,fi,x);

    hold on
    plot (x,f(x),'k') %real
    plot (x,p)
    hold off
    a=c; %avanzar
end

end
%% 2 
% hecho en la hoja

%% 3a

f= @(x) atan(100.*(x-0.873));

a = -3;
b = 2;

Ne = 28;          % Número de elementos 
m=1;
x = a:0.01:b;

figure ('Name', 'Ne=28 m=1')
p=elementos_finitos_1D(f,a,b,Ne,m);
error1 = abs(max(f(x)) - max(p)); % para calcular el error elegimos el punto máximo de las funciones en el
                                    % intervalo dado (en otros puntos seguramente el error cambie)

Ne = 14;         
m=2;
x = a:0.01:b
figure ('Name', 'Ne=14 m=2')
p=elementos_finitos_1D(f,a,b,Ne,m);
error2 = abs(max(f(x)) - max(p));


Ne = 7;         
m=4;
x = a:0.01:b
figure ('Name', 'Ne=7 m=4')
p=elementos_finitos_1D(f,a,b,Ne,m);
error3 = abs(max(f(x)) - max(p));

Ne = 1;         
m=28;
x = a:0.01:b
figure ('Name', 'Ne=1 m=28')
p=elementos_finitos_1D(f,a,b,Ne,m);
error4= abs(max(f(x)) - max(p));

error1
error2
error3
error4


%vemos que el error más bajo se da con Ne = 14 y m=2  error===1.026698484896116e-06



%% 3b

format long

%delta = 1e-2

a = -3;
b = 2;
long = b-a;
delta = 1e-2;
h1 = zeros(4,1); 
Ne1 = zeros(4,1); 
Nn1 = zeros(4,1); 

x=-1:0.01:1;

xi=linspace(-1,1,2);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC1=max(abs(w));

xi=linspace(-1,1,3);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC2=max(abs(w));

xi=linspace(-1,1,5);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC4=max(abs(w));

xi=linspace(-1,1,29);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC28=max(abs(w));

cot_der = [10^(3), 2*10^(6), 2.4*10^(11), 3.1*10^(87)]; 
O = [2, 3, 5, 29]; % m + 1
C = [maximoC1, maximoC2, maximoC4, maximoC28];

for i = 1:4
 h1(i) = ((delta*factorial(O(i))*2^(O(i))/(cot_der(i)*C(i))))^(1/O(i));
 Ne1 (i) = long/h1(i);
 Nn1 (i) = (O(i)-1)*Ne1(i)+1;
end
h1
Ne1
Nn1

% delta = 1e-6;
format long
a = -3;
b = 2;
long = b-a;
delta = 1e-6;
h2 = zeros(4,1); 
Ne2 = zeros(4,1); 
Nn2 = zeros(4,1); 

x=-1:0.01:1;

xi=linspace(-1,1,2);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC1=max(abs(w));

xi=linspace(-1,1,3);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC2=max(abs(w));

xi=linspace(-1,1,5);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC4=max(abs(w));

xi=linspace(-1,1,29);
w=1;
for i=1 : length(xi)
    w = w.*(x-xi(i));
end
maximoC28=max(abs(w));

cot_der = [10^(3), 2*10^(6), 2.4*10^(11), 3.1*10^(87)]; 
O = [2, 3, 5, 29]; % m + 1
C = [maximoC1, maximoC2, maximoC4, maximoC28];

for i = 1:4
 h2(i) = ((delta*factorial(O(i))*2^(O(i))/(cot_der(i)*C(i))))^(1/O(i));
 Ne2 (i) = long/h2(i);
 Nn2 (i) = (O(i)-1)*Ne2(i)+1;
end
h2
Ne2
Nn2




%% EJERCICIO 2


%creación del mallado
[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);

%matrices donde se van a almacenar los puntos transformados
   %almacenaje de las coordenadas x
coordx=zeros(11,11);
coordy=zeros(11,11);

%nodos del cuadrilátero distorsionado
X1=[0;0];
X2=[0.6;0];
X3=[0.5;0.5];
X4=[0;0.3];
X=[X1 X2 X3 X4]';  %trasponer para que coincidan dimensiones


%para realizar la transformación sabemos
%(x,y)=Fi(xg,yg)=(sum(xj*phij),sum(yj*phij) para todo m asi que aplicamos esta fórmula

for i=1:11
    for j=1:11  %evaluamos las funciones base en el nodo (j,i) del mallado de referencia
        phi1=0.25*(xg(j,i)-1)*(yg(j,i)-1);
        phi2=-0.25*(xg(j,i)+1)*(yg(j,i)-1);
        phi3=0.25*(xg(j,i)+1)*(yg(j,i)+1);
        phi4=-0.25*(xg(j,i)-1)*(yg(j,i)+1);

        phi=[phi1;phi2;phi3;phi4];
        for k=1:4  %almacenaje en las matrices siguiendo la fórmula de (x,y)=Fi(xg,yg)
            coordx(j,i)=coordx(j,i)+X(k,1)*phi(k);
            coordy(j,i)=coordy(j,i)+X(k,2)*phi(k);
        end
    end
end


%Apartado 2)

f=@(x,y) sin(pi*x)*sin(pi*y);

%creación del mallado
[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);

%matrices donde se van a almacenar los puntos transformados
   %almacenaje de las coordenadas x
coordx=zeros(11,11);
coordy=zeros(11,11);

%nodos del cuadrilátero distorsionado
X1=[0;0];
X2=[0.6;0];
X3=[0.5;0.5];
X4=[0;0.3];
X=[X1 X2 X3 X4]';  %trasponer para que coincidan dimensiones

f1=f(X1(1),X1(2));
f2=f(X2(1),X2(2));
f3=f(X3(1),X3(2));
f4=f(X4(1),X4(2));


%para realizar la transformación sabemos
%(x,y)=Fi(xg,yg)=(sum(xj*phij),sum(yj*phij) asi que aplicamos esta fórmula
%también sabemos que fh=f1+phi1+f2*phi2+f3*phi3+f4*phi4
for i=1:11
    for j=1:11  %evaluamos las funciones base en el nodo (j,i) del mallado de referencia
        phi1=0.25*(xg(j,i)-1)*(yg(j,i)-1);
        phi2=-0.25*(xg(j,i)+1)*(yg(j,i)-1);
        phi3=0.25*(xg(j,i)+1)*(yg(j,i)+1);
        phi4=-0.25*(xg(j,i)-1)*(yg(j,i)+1);

        phi=[phi1;phi2;phi3;phi4];
        for k=1:4  %almacenaje en las matrices siguiendo la fórmula de (x,y)=Fi(xg,yg)
            coordx(j,i)=coordx(j,i)+X(k,1)*phi(k);
            coordy(j,i)=coordy(j,i)+X(k,2)*phi(k);
        end
        fh(j,i)=f1*phi1+f2*phi2+f3*phi3+f4*phi4;  
    end
end

%dibujar superficie
surf(coordx,coordy,fh)



%Apartado 3)
%Hay que repetirlo cambiando las coordenadas del cuadrilatero, es decir
%cambiando X1,X2,X3,X4
figure 
surf(coordx,coordy,fh)
hold on

%CUADRILÁTERO 2
f=@(x,y) sin(pi*x)*sin(pi*y);

%creación del mallado
[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);

%matrices donde se van a almacenar los puntos transformados
   %almacenaje de las coordenadas x
coordx=zeros(11,11);
coordy=zeros(11,11);

%nodos del cuadrilátero distorsionado
X1=[0.6;0];
X2=[1;0];
X3=[1;0.4];
X4=[0.5;0.5];
X=[X1 X2 X3 X4]';  %trasponer para que coincidan dimensiones

f1=f(X1(1),X1(2));
f2=f(X2(1),X2(2));
f3=f(X3(1),X3(2));
f4=f(X4(1),X4(2));


%para realizar la transformación sabemos
%(x,y)=Fi(xg,yg)=(sum(xj*phij),sum(yj*phij) asi que aplicamos esta fórmula
%también sabemos que fh=f1+phi1+f2*phi2+f3*phi3+f4*phi4
for i=1:11
    for j=1:11  %evaluamos las funciones base en el nodo (j,i) del mallado de referencia
        phi1=0.25*(xg(j,i)-1)*(yg(j,i)-1);
        phi2=-0.25*(xg(j,i)+1)*(yg(j,i)-1);
        phi3=0.25*(xg(j,i)+1)*(yg(j,i)+1);
        phi4=-0.25*(xg(j,i)-1)*(yg(j,i)+1);

        phi=[phi1;phi2;phi3;phi4];
        for k=1:4  %almacenaje en las matrices siguiendo la fórmula de (x,y)=Fi(xg,yg)
            coordx(j,i)=coordx(j,i)+X(k,1)*phi(k);
            coordy(j,i)=coordy(j,i)+X(k,2)*phi(k);
        end
        fh2(j,i)=f1*phi1+f2*phi2+f3*phi3+f4*phi4;  
    end
end
surf(coordx,coordy,fh2)
%CUADRILÁTERO 3
f=@(x,y) sin(pi*x)*sin(pi*y);

%creación del mallado
[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);

%matrices donde se van a almacenar los puntos transformados
   %almacenaje de las coordenadas x
coordx=zeros(11,11);
coordy=zeros(11,11);

%nodos del cuadrilátero distorsionado
X1=[0.5;0.5];
X2=[1;0.4];
X3=[1;1];
X4=[0.5;1];
X=[X1 X2 X3 X4]';  %trasponer para que coincidan dimensiones

f1=f(X1(1),X1(2));
f2=f(X2(1),X2(2));
f3=f(X3(1),X3(2));
f4=f(X4(1),X4(2));


%para realizar la transformación sabemos
%(x,y)=Fi(xg,yg)=(sum(xj*phij),sum(yj*phij) asi que aplicamos esta fórmula
%también sabemos que fh=f1+phi1+f2*phi2+f3*phi3+f4*phi4
for i=1:11
    for j=1:11  %evaluamos las funciones base en el nodo (j,i) del mallado de referencia
        phi1=0.25*(xg(j,i)-1)*(yg(j,i)-1);
        phi2=-0.25*(xg(j,i)+1)*(yg(j,i)-1);
        phi3=0.25*(xg(j,i)+1)*(yg(j,i)+1);
        phi4=-0.25*(xg(j,i)-1)*(yg(j,i)+1);

        phi=[phi1;phi2;phi3;phi4];
        for k=1:4  %almacenaje en las matrices siguiendo la fórmula de (x,y)=Fi(xg,yg)
            coordx(j,i)=coordx(j,i)+X(k,1)*phi(k);
            coordy(j,i)=coordy(j,i)+X(k,2)*phi(k);
        end
        fh3(j,i)=f1*phi1+f2*phi2+f3*phi3+f4*phi4;  
    end
end
surf(coordx,coordy,fh3)
%CUADRILÁTERO 4
f=@(x,y) sin(pi*x)*sin(pi*y);

%creación del mallado
[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);

%matrices donde se van a almacenar los puntos transformados
   %almacenaje de las coordenadas x
coordx=zeros(11,11);
coordy=zeros(11,11);

%nodos del cuadrilátero distorsionado
X1=[0;0.3];
X2=[0.5;0.5];
X3=[0.5;1];
X4=[0;1];
X=[X1 X2 X3 X4]';  %trasponer para que coincidan dimensiones

f1=f(X1(1),X1(2));
f2=f(X2(1),X2(2));
f3=f(X3(1),X3(2));
f4=f(X4(1),X4(2));


%para realizar la transformación sabemos
%(x,y)=Fi(xg,yg)=(sum(xj*phij),sum(yj*phij) asi que aplicamos esta fórmula
%también sabemos que fh=f1+phi1+f2*phi2+f3*phi3+f4*phi4
for i=1:11
    for j=1:11  %evaluamos las funciones base en el nodo (j,i) del mallado de referencia
        phi1=0.25*(xg(j,i)-1)*(yg(j,i)-1);
        phi2=-0.25*(xg(j,i)+1)*(yg(j,i)-1);
        phi3=0.25*(xg(j,i)+1)*(yg(j,i)+1);
        phi4=-0.25*(xg(j,i)-1)*(yg(j,i)+1);

        phi=[phi1;phi2;phi3;phi4];
        for k=1:4  %almacenaje en las matrices siguiendo la fórmula de (x,y)=Fi(xg,yg)
            coordx(j,i)=coordx(j,i)+X(k,1)*phi(k);
            coordy(j,i)=coordy(j,i)+X(k,2)*phi(k);
        end
        fh4(j,i)=f1*phi1+f2*phi2+f3*phi3+f4*phi4;  
    end
end
surf(coordx,coordy,fh4)
hold off











