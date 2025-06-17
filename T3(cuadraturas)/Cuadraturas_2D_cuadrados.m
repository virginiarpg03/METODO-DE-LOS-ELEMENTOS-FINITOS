%% CUADRATURAS EN 2D (para cuadrados)

clear all
close all
clc 

% Definir la función f(x, y)
f = @(x, y) x.^3 + y.^3;  % Ejemplo de función en 2D

%definir los vértices de tal forma que se correspondan las coordenadas x e y 
Vx = [];  
Vy = [];  


%% Regla de los trapecios en 2D
Int = zeros(--, 1);  % Vector para almacenar las integrales 

for i = 1:16
    % coordenadas
    x1 = Vx(i, 1); x2 = Vx(i, 2);  hx = (x2 - x1);
    y1 = Vy(i, 2); y2 = Vy(i, 3);  hy = (y2 - y1);
    
    Int(i) = hx * hy / 4 * (f(x1, y1) + f(x1, y2) + f(x2, y1) + f(x2, y2));  % hx * hy es el área del cuadrado
end

Q_trapecios = sum(Int);

%% Regla de Simpson en 2D
Int = zeros( , 1);  

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

%para cuadriláteros deformados necesitamos el jacobiano para el área (ej 8)

%% Cuadratura de Gauss-Legendre 
% Definir puntos y pesos de Gauss-Legendre en [-1, 1]

xi = [ ];
wi = [ ];

Int = zeros(16, 1);  % Vector para almacenar las integrales de cada subdominio

for i = 1:16
    % Extraer las coordenadas de los vértices
    x1 = Vx(i, 1); x2 = Vx(i, 2);  hx = (x2 - x1);
    y1 = Vy(i, 2); y2 = Vy(i, 3);  hy = (y2 - y1);
    

 Intx=hx /2*sum(f1((xi+1).*(x2-x1)+x1).*wi);
 Inty=hy/2*sum(f2((xi+1).*(y2-y1)+y1).*wi);
 Int(i)=Intx*Inty;
end

Q_gauss = sum(Int);

%realmente si mis cuadrados no están entre 0 y 1 me suicido