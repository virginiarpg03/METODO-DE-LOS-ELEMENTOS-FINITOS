%% CUADRATURAS EN 2D (para cuadrados)

clear all
close all
clc 

f = @(x,y) sin(pi*x).*sin(pi*y);

% el malladao lo hago por cuadrados y meto los vértices a mano :(

Vx= [0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75]; 

Vy= [0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1]; 

%% Regla de los trapecios en 2D
Int = zeros(16, 1);  % Vector para almacenar las integrales 

for i = 1:16
    % coordenadas
    x1 = Vx(i, 1); x2 = Vx(i, 2);  hx = (x2 - x1);
    y1 = Vy(i, 2); y2 = Vy(i, 3);  hy = (y2 - y1);
    
    Int(i) = hx * hy / 4 * (f(x1, y1) + f(x1, y2) + f(x2, y1) + f(x2, y2));  % hx * hy es el área del cuadrado
end

Q_trapecios = sum(Int);

%% Regla de Simpson en 2D          algo mal?¿
Int = zeros(16 , 1);  

for i = 1:16
    % coordenadas de los vértices
    x1 = Vx(i, 1); x2 = Vx(i, 2);  hx = (x2 - x1);
    y1 = Vy(i, 2); y2 = Vy(i, 3);  hy = (y2 - y1);
    
    % Puntos medios en x e y
    xm = (x1 + x2) / 2;
    ym = (y1 + y2) / 2;
    
    % Aplicar la regla de Simpson en 2D
    Int(i) = hx * hy / 9 * (f(x1, y1) + 4 * f(xm, y1) + f(x2, y1) + ...
             4 * f(x1, ym) + 16 * f(xm, ym) + 4 * f(x2, ym) +  f(x1, y2) + 4 * f(xm, y2) + f(x2, y2));
end

Q_simpson = sum(Int);

%% Cuadratura de Gauss-Legendre 
% Definir puntos y pesos de Gauss-Legendre en [-1, 1]

xi=[-0.8611363115940526 -0.3399810435848563 0.3399810435848563 0.8611363115940526];
wi=[0.3478548451374538 0.6521451548625461 0.6521451548625461 0.3478548451374538];

Int = zeros(16, 1);  % Vector para almacenar las integrales de cada subdominio
f1=@(x) sin(pi*x);
f2=@(y) sin(pi*y);
for i = 1:16
    % Extraer las coordenadas de los vértices
    x1 = Vx(i, 1); x2 = Vx(i, 2);  hx = (x2 - x1);
    y1 = Vy(i, 2); y2 = Vy(i, 3);  hy = (y2 - y1);
    

 Intx=hx /2*sum(f1((xi+1).*(x2-x1)+x1).*wi);
 Inty=hy/2*sum(f2((xi+1).*(y2-y1)+y1).*wi);
 Int(i)=Intx*Inty;
end

Q_gauss = sum(Int);

%% Resultados
disp(['Integral con la regla de los trapecios: ', num2str(Q_trapecios)]);
disp(['Integral con la regla de Simpson: ', num2str(Q_simpson)]);
disp(['Integral con la cuadratura de Gauss-Legendre: ', num2str(Q_gauss)]);

%%


Int=zeros(16,1); %Vector de ceros para la integral en cada elemento
for i=1:16
 x1=Vx(i,1); x2=Vx(i,2); y1=Vy(i,2); y2=Vy(i,3);
 Int(i)=(x2-x1)*(y2-y1)/4*(f(x1,y1)+f(x1,y2)+f(x2,y1)+f(x2,y2));
end
I_calc=sum(Int)