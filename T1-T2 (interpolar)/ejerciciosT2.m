%% EJERCICIOS HOJA 2



%% ejercicio 6


% a) % polinomios lineales: 

% el triángulo que contiene al punto que buscamos es el de coordenadas:
% P1[0.25;0.5] P2 [0.5;0.75]     P3[0.25;0.75]

f=@(x,y) sin(pi*x).*sin(pi*y);    %funcion a interpolar
X=[0.32;0.74]; %punto a interpolar

%coordenadas triángulo que contiene al punto
X1 = [0.25;0.5];
X2 = [0.5;0.75];     
X3 = [0.25;0.75];

% buscar xgorro, ygorro
bi=X1;
Ai= [X2-X1 , X3-X1];
Xg= Ai\(X-bi);

xg=Xg(1);
yg=Xg(2);

% queremos hallar phi(X). Necesitamos la phi gorro ya que ==> phi(X)=phigorro(Xg) 

f1= f(X1(1) , X1(2));
f2= f(X2(1) , X2(2));
f3= f(X3(1) , X3(2));

phi1= 1 -xg -yg;
phi2= xg;
phi3= yg;

fh= f1*phi1 + f2*phi2 + f3*phi3

% CON LA FUNCION
fh_lineal = elementos_triangulos_lineal(f,X1,X2,X3,X)


%b) polinomios cuadráticos

%ahora tambien los puntos medios

fh_cuadratico = elementos_triangulos_cuadratico(f,X1,X2,X3,X)

%% ejercicio 7

f=@(x,y) exp(-((x-1).^2 + y.^2)./0.4^2);

% he utilizado pdtool para hacer una circunc¡ferencia centrada en (1,0) y
% con radio 1; esportando :
% P == coordenadas de los nodos
% t== orden de los nodos

t=t(1:3,:);      %la fila 4 no
X=[1;0];        %punto a interpolar
k=tsearchn(p',t', [1 0]);            % nº de la columna de la matriz t que contiene los vértices del 
                                    % triángulo que contiene al punto X (las coordenadas guardadas en p)


% vertices del triangulo:                                    
X1 = p(:,t(1,k));
X2 = p(:,t(2,k));
X3 = p(:,t(3,k));

fh_lineal = elementos_triangulos_lineal(f,X1,X2,X3,X)

fh_cuadra = elementos_triangulos_cuadratico(f,X1,X2,X3,X)




%% ejercicio 8


f=@(x,y) sin(pi*x).*sin(pi*y);
X=[0.32 ; 0.74];

%VÉRTICES CUADRADO
V1=[0.25; 0.5];
V2=[0.5; 0.5];
V3=[0.5; 0.75];
V4=[0.25; 0.75];

fh_lineal = elementos_cuadrados(f,V1,V2,V3,V4,X,1)
fh_cuadratico = elementos_cuadrados(f,V1,V2,V3,V4,X,2)


%% ejercicio 9


format long

a = -3;
b = 2;
long = b - a;
delta = 1e-6;

h = zeros(4,1); %vector de pasos
Ne = zeros(4,1); %vector de numero de elementos
Nn = zeros(4,1); %vector de numero de nodos

c = [10^3, 2*10^6, 2.4*10^11, 3.1*10^87]; %cotas de derivadas
B = [2, 3, 5, 29]; %grado del polinomio + 1

for i = 1:4
    h(i) = ((delta/(c(i))))^(1/B(i));
    Ne(i) = long/h(i);
    Nn(i) = (B(i)-1)*Ne(i)+1;
end

h