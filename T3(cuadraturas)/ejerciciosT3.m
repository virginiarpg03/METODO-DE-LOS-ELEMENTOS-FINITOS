%% EJERCICIOS T3


%% EJ 1

% a)
 %resolvemos el sistema de ecuaciones (cuadratura exacta de 1 , x , x^2)

A=[1 1 1; 1/4 1/2 3/4; 1/16 1/4 9/16];
b=[1; 1/2; 1/3];
wi= A\b;

% b)
f=@(x) x.^2/((x.^2+12).^(1/2));


xi=[1/4 1/2 3/4];
Qf=sum(w*f(xi));

%% EJ 2

% Sea Dh una partición equiespaciada del intervalo [1, 2]. Se pide calcular el valor aproximado de la integral
% mediante las siguientes cuadraturas:
f=@(x) 1./x;
a= 1;
b=2;
h= [0.5 0.25 0.175 0.0625];

Qf_trapecios=zeros(length(h),1);
Qf_simpsom = Qf_trapecios;
Qf_legendre4 = Qf_trapecios;
Qf_legendre8 = Qf_trapecios;

% a)  regla de los trapecios compuesta.

for k=1:length(h)
    n= (b-a)/h(k);
    sumf= 0;

    for i= 1:n-1
        x(i)= a+h(k)*i;
        sumf= sumf + f(x(i));
    end
    
    Qf_trapecios(k)= h(k)*0.5*(f(a) + 2*sumf + f(b));

end

Qf_trapecios

% b) simpsom compuesta

for k=1:length(h)
    n= (b-a)/h(k);
    s1=0;
    s2=0;

    for i=1:n-1
        x=a +i*h(k);
    
        if rem(i,2)==0
           s2 = s2 + f(x);
        else
           s1 = s1 + f(x);
        end
    end
    
    Qf_simpsom(k) = h(k)/3 * (f(a) + 4*s1 + 2*s2 + f(b));

end

Qf_simpsom


% c)  regla del punto medio compuesta. no lo hago

% d) cuadratura compuesta de Gauss-Legendre con n = 4 nodos.

xi = [-0.8611363115940526 -0.3399810435848563 0.3399810435848563 0.8611363115940526];
wi = [0.3478548451374538 0.6521451548625461 0.6521451548625461 0.3478548451374538];

for k = 1:length(h)
    n = round((b - a) / h(k));  % Aseguramos que n sea un número entero
    zi = linspace(a, b, n+1); 
    Qf_legendre4(k) = 0;  

    for i = 1:n
        ai = zi(i);
        bi = zi(i+1);
        
        % Definimos la función g en el subintervalo [a_i, b_i] utilizando la transformación lineal
        g = @(x) 0.5*(bi - ai) * f(0.5*(bi - ai)*(x + 1) + ai);
        
        % Sumar el valor de la cuadratura de Gauss-Legendre en el subintervalo [a_i, b_i]
        Qf_legendre4(k) = Qf_legendre4(k) + sum(wi .* g(xi));
    end
end

Qf_legendre4


% e) cuadratura compuesta de Gauss-Legendre con n = 8 nodos.

xi = [-0.9602898564975363 -0.7966664774136267 -0.5255324099163290 -0.1834346424956498 0.1834346424956498 0.5255324099163290 0.7966664774136267 0.9602898564975363];
wi = [0.1012285362903763 0.2223810344533745 0.3137066458778873 0.3626837833783620 0.3626837833783620 0.3137066458778873 0.2223810344533745 0.1012285362903763];

for k = 1:length(h)
    n = round((b - a) / h(k));  % Aseguramos que n sea un número entero
    zi = linspace(a, b, n+1); 
    Qf_legendre8(k) = 0;  

    for i = 1:n
        ai = zi(i);
        bi = zi(i+1);
        
        % Definimos la función g en el subintervalo [a_i, b_i] utilizando la transformación lineal
        g = @(x) 0.5*(bi - ai) * f(0.5*(bi - ai)*(x + 1) + ai);
        
        % Sumar el valor de la cuadratura de Gauss-Legendre en el subintervalo [a_i, b_i]
        Qf_legendre8(k) = Qf_legendre8(k) + sum(wi .* g(xi));
    end
end

Qf_legendre8


%los errores los va acalcular su madre la coja que pesadilla




%% EJ 6 cuadraditos

%creo el mallado y defino la función

f = @(x,y) sin(pi*x).*sin(pi*y);

% el malladao lo hago por cuadrados y meto los vértices a mano :(

   %     4           3



   %     1           2

Vx= [0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75;...
     0 .25 .25 0; .25 .5 .5 .25; .5 .75 .75 .5; .75 1 1 .75]; 

Vy= [0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1;...
     0 0 .25 .25; .25 .25 .5 .5; .5 .5 .75 .75; .75 .75 1 1]; 

%son matrices de 16x4 (16 cuadrados x 4 vértices cada uno)

% a)
Int = zeros(16,1);  %valor de la integral en cada elemento

for i=1:16
        % coordenadas del cuadrado que integro

        x1 = Vx(i,1); x2 = Vx(i,2);  hx=(x2 - x1) ;
        y1 = Vy(i,2); y2 = Vy(i,3);  hy=(y2 - y1) ;
        
        % Regla de los trapecios en 2D 
        Int(i) = hx*hy  / 4 * (f(x1, y1) + f(x1, y2) + f(x2, y1) + f(x2, y2) ); %hx*hy es el área de cada cuadrado
       
        
 end

 Q_trapecio= sum(Int)


 %% EJ 8

f = @(x,y) sin(pi*x).*sin(pi*y);  

% Vértices de los cuadriláteros (cada columna representa un cuadrilátero)
nx = [0 0.6 0.5 0; 0.6 1 1 0.5; 0.5 1 1 0.5; 0 0.5 0.5 0];
ny = [0 0 0.5 0.3; 0 0 0.4 0.5; 0.3 0.5 1 1; 0.5 0.4 1 1];

% cuadrilátero de referencia [-1, 1] x [-1, 1]
Xg = [-1; 1; 1; -1; 0; 1; 0; -1; 0];  
Yg = [-1; -1; 1; 1; -1; 0; 1; 0; 0]; 
w = [1; 1; 1; 1; 4; 4; 4; 4; 16]; 

% Funciones 
phi1 = @(x, y) 1/4 * (1 - x) .* (1 - y);  
phi2 = @(x, y) 1/4 * (1 + x) .* (1 - y); 
phi3 = @(x, y) 1/4 * (1 + x) .* (1 + y);
phi4 = @(x, y) 1/4 * (1 - x) .* (1 + y);  

% Inicializar la variable para acumular la integral total
q = 0;

% Bucle sobre los cuadriláteros
for i = 1:4
    % Extraer las coordenadas de los vértices del cuadrilátero i
    xn = nx(:, i);
    yn = ny(:, i);

    % Definir las funciones de coordenadas del cuadrilátero transformado
    x_real = @(x, y) xn(1)*phi1(x, y) + xn(2)*phi2(x, y) + xn(3)*phi3(x, y) + xn(4)*phi4(x, y);
    y_real = @(x, y) yn(1)*phi1(x, y) + yn(2)*phi2(x, y) + yn(3)*phi3(x, y) + yn(4)*phi4(x, y);

    % Bucle sobre los puntos de cuadratura (9 puntos)
    for j = 1:9
        % Calcular las coordenadas (x, y) transformadas en el cuadrilátero real
        x = x_real(Xg(j), Yg(j));
        y = y_real(Xg(j), Yg(j));

        % Cálculo del Jacobiano
        dxdx = 1/4 * ((1 - Yg(j)) * (xn(2) - xn(1)) + (1 + Yg(j)) * (xn(3) - xn(4)));
        dxdy = 1/4 * ((1 - Xg(j)) * (xn(4) - xn(1)) + (1 + Xg(j)) * (xn(3) - xn(2)));
        dydx = 1/4 * ((1 - Yg(j)) * (yn(2) - yn(1)) + (1 + Yg(j)) * (yn(3) - yn(4)));
        dydy = 1/4 * ((1 - Xg(j)) * (yn(4) - yn(1)) + (1 + Xg(j)) * (yn(3) - yn(2)));

        % Determinante del Jacobiano (cambio de área)
        JF = dxdx * dydy - dxdy * dydx;

        % Corregir el uso del Jacobiano absoluto
        q = q + w(j) * f(x, y) * abs(JF);  % El Jacobiano debe ser positivo
    end
end

% Resultado final de la integral
qfinal = q / 9;  % Dividimos por 9 por la regla de Simpson

% Mostrar el resultado
disp(['Integral final: ', num2str(qfinal)]);


%% EJ 9


load('malla_ej9.mat')

f=@(x,y) exp(-((x-1).^2 +y.^2)/0.4^2);

xi=p(1,:);    %valor de la x de cada nodo
yi=p(2,:);    %valor de la y de cada nodo

elem=t(1:3,:)';    %guarda los nodos de cada triángulo

trisurf(elem,xi,yi,f(xi,yi)) %dibujo


Ne=length(elem);    %nos da elnúmero de triángulos que ha generado el mallado

Xg=[0 1 0;0 0 1];    %los valores de (xg,yg) que se usan para la cuadratura
wi=[1/6 1/6 1/6];    %pesos que se usan en la cuadratura

Q=0;

for i=1:Ne    %las coordenadas de cada nodo de cada triángulo
    X1=p(:,elem(i,1));
    X2=p(:,elem(i,2));
    X3=p(:,elem(i,3));

%ahora aplicamos la transformación afín para triángulos 
bi=X1;
Ai=[(X2-X1) (X3-X1)];

    for j=1:3
        X=bi+Ai*Xg(:,j); %creamos los tres vértices transformados 
        Q=Q+wi(j)*abs(det(Ai))*f(X(1), X(2));
     end 
end 
Q


%%
