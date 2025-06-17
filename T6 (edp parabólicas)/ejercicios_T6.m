%% EJERCICIOS T6


%% REGIÓN ESTABILIDAD    dt < 2/(k * Dmax)

% Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A = M;
A0 = A;
A0([1 end],:) = 0;
A0(:,[1 end]) = 0;
A0(1,1) = 1;
A0(end,end) = 1;

R0 = R;             
R0([1 end],:) = 0;
R0(:,[1 end]) =0;
R0(1,1) = 1;
R0(end,end) = 1;

% autovalores (D) y autovectores (V) del problema
[V D]=eig(full(R0),full(M));



%% EJERCICIO 1 : EULER EXPLÍCITO

clear all

a= 0;
k=0.1;

% intervalo espacio:
c = 0;
d = 1;

Ne = 100;           % número elementos
h = (d-c)/Ne;
xi = c:h:d;         % nodos

% intervalo tiempo:
t0= 0;
tf= 100;
dt= 0.00016;                % para que el euler explícito sea estable: dt<h^2/2kpi
Nt= (tf-t0)/dt;

% Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A = M;
A0 = A;
A0([1 end],:) = 0;
A0(:,[1 end]) = 0;
A0(1,1) = 1;
A0(end,end) = 1;

uhn = sin(pi*xi').^100;
plot(xi,uhn)
title(0)
pause(0.0005)

% resolver el pvi con euler explícito 
for n = 1:Nt
    vect_b = ((1-a*dt)*M-k*dt*R)*uhn;
    vect_b0 = vect_b;
    vect_b0([1 end])=0;

    uhn = A0\vect_b0; 

    plot(xi,uhn)
    title(n*dt)
    axis([0 1 -.2 1.2])

    pause(0.00005)
    %Ctrol+C para parar la grafica
end

%% EJERCICIO 2 : EULER IMPLÍCITO

clear all

a= 0;
k=0.1;

% intervalo espacio:
c = 0;
d = 1;

Ne = 10000;           % número elementos
h = (d-c)/Ne;
xi = c:h:d;         % nodos

% intervalo tiempo:
t0= 0;
tf= 100;
dt= 0.01;                % vemos como el euler implícito es incondicionalmente estable
Nt= (tf-t0)/dt;

% Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A = (M + R*k*dt);
A0 = A;
A0([1 end],:) = 0;
A0(:,[1 end]) = 0;
A0(1,1) = 1;
A0(end,end) = 1;

uhn = sin(pi*xi').^100;
plot(xi,uhn)
title(0)
pause(3)

% resolver el pvi con euler implícito 
for n = 1:Nt
    vect_b = M*uhn;
    vect_b0 = vect_b
    vect_b0([1 end])=0;

    uhn = A0\vect_b0; 

    plot(xi,uhn)
    title(n*dt)
    axis([0 1 -.2 1.2])

    pause(0.0005)
    %Ctrol+C para parar la grafica
end


%% EJERCICIO 3 : TRAPECIOS (CRANK-NICOLSON)

clear all

a= 0;
k=0.1;

% intervalo espacio:
c = 0;
d = 1;

Ne = 100;           % número elementos
h = (d-c)/Ne;
xi = c:h:d;         % nodos

% intervalo tiempo:
t0= 0;
tf= 100;
dt= 0.001;                % vemos como es A-estable
Nt= (tf-t0)/dt;

% Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A = (M + R*k*dt*0.5);
A0 = A;
A0([1 end],:) = 0;
A0(:,[1 end]) = 0;
A0(1,1) = 1;
A0(end,end) = 1;

uhn = sin(pi*xi').^100;
plot(xi,uhn)
title(0)
pause(3)

% resolver el pvi con trapecios
for n = 1:Nt
    vect_b = (M - k*dt*0.5*R)*uhn;
    vect_b0 = vect_b
    vect_b0([1 end])=0;

    uhn = A0\vect_b0; 

    plot(xi,uhn)
    title(n*dt)
    axis([0 1 -.2 1.2])

    pause(0.0005)
    %Ctrol+C para parar la grafica
end

%% EJERCIIO 4 : edp parabólica con trapecios y condiciones no homo

clear;

% Definir la función 
f = @(x) 50 * exp(-x.^2 / 0.01);
a = 0; 
k = 0.01;

% Intervalo espacio
c = -1; 
d = 4;
Ne = 1000;               % Número de elementos
h = (d - c) / Ne;        % Tamaño del paso 
xi = c:h:d;              % Nodos del mallado 


% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.016;           
Nt = (tf - t0) / dt;    

 % Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A = M + 0.5*dt*k*R;
A0 = A;
A0([1, end], :) = 0;     
A0(:, [1, end]) = 0;    
A0(1, 1) = 1;           
A0(end, end) = 1;       

% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(1) = 5;               
gi(end) = 30;      

% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 20;
whn = uhn;
fi = f(xi)';

plot(xi', uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt

    vect_b = (M - 0.5*k*dt*R)*uhn + dt*M*(fi) - A*gi;
    vect_b0 = vect_b;
    vect_b0(1) = 0;          
    vect_b0(end) = 0;       

    % Resolver el sistema lineal
    whn = A0 \ vect_b0;
    uhn = whn + gi;
    
    plot(xi,uhn)
    title(n*dt)

     pause
end



%% EJERCICIO 5 : condiciones de contorno no homo y f(x,t)

clear;

% Definir la función 
f = @(x,t) 50 * exp((-x.^2 / 0.01) - (t-20.^2)) ;
 a = 0; 
k = 0.01;

% Intervalo espacio
c = -1; 
d = 4;
Ne = 100;               % Número de elementos
h = (d - c) / Ne;        % Tamaño del paso 
xi = c:h:d;              % Nodos del mallado 


% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.05;           
Nt = (tf - t0) / dt;    

% Matrices de masas y rigidez las cuales no varían
Calcular_Matriz_Masa_Rigidez_Lineal1D

A= (1+a*dt*0.5)*M + dt*k*0.5*R;

A0 = A; 
A0(end,:) = 0; 
A0(:,end) = 0; 
% A0(1,1) = 1;       % NO CONOCEMOS LA CONDICION DE CONTORNO EN d
A0(end,end) = 1;  


% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(end) = 30;      

% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 20;
whn = uhn;
plot(xi', uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt
    fi0= f(xi', (n-1)*dt);
    fi1= f(xi', n*dt);
    vect_b = ((1-a*dt*0.5)*M-k*dt*0.5*R)*uhn + dt*0.5*M*(fi0 + fi1)-A*gi;
    %Condiciones de contorno vect_b -> vect_b0
    vect_b0=vect_b;
    vect_b0(1)=vect_b0(1)+6*dt*k;
    vect_b0(end)=0;

    whn = A0\vect_b0;
    uhn= whn+gi;

    plot(xi,uhn,'.-')
    title(n*dt)
    pause(0.05)
end


%% EJERCICIO 6

clear;

% Definir la función 
f = @(x,y) 0.*x + 0.*y;
a = 0; 
k = 0.01;

    % MALLADO
load("malla_ej6.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = unique([e(1,:) e(2,:)]); %no repetidos
ind_1 = find(yi==1);

% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.2;           
Nt = (tf - t0) / dt;    

% Matrices de masas y rigidez las cuales no varían
[R M] = assema(p,t,1,1,0);

A = M + 0.5*dt*k*R;
A0 = A;
A0(fron_d, :) = 0;     
A0(:, fron_d) = 0;    
for i = fron_d
    A0(i,i) = 1;
end
      

% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(ind_1) = 50;               

% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 0*yi' + 30;
whn = uhn;
fi = f(xi,yi)';

trisurf(elem,xi, yi, uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt

    vect_b = (M - 0.5*k*dt*R)*uhn + dt*M*(fi) - A*gi;
    vect_b0 = vect_b;
    vect_b0(fron_d) = 0;          

    % Resolver el sistema lineal
    whn = A0 \ vect_b0;
    uhn = whn + gi;
    
trisurf(elem,xi,yi,uhn)
    title(n*dt)

     pause
end

%% EJERCICIO 7 : la frontera dirichlet cambia en un t

clear;

% Definir la función 
f = @(x,y) 0.*x + 0.*y;
a = 0; 
k = 0.01;

    % MALLADO
load("malla_ej6.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_t = unique([e(1,:) e(2,:)]); %no repetidos
fron1 = find(yi==1);
fron2 = find(yi==0 | xi==0 | xi==1);

% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.2;           
Nt = (tf - t0) / dt;    

% Matrices de masas y rigidez 
[R M] = assema(p,t,1,1,0);

A = M + 0.5*dt*k*R;

% Entre 0 y 10 la frontera dirichlet será toda la frontera
A01 = A;
A01(fron_t, :) = 0;     
A01(:, fron_t) = 0;    
for i = fron_t
    A01(i,i) = 1;
end

% Entre 10 y 100 la frontera drichlet solo será la frontera 1
A02 = A;
A02(fron1, :) = 0;     
A02(:, fron1) = 0;    
for i = fron1
    A02(i,i) = 1;
end
      

% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(fron1) = 50;               

% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 0*yi' + 30;
whn = uhn;
fi = f(xi,yi)';

trisurf(elem,xi, yi, uhn);
pause;

% resolver problemas elípticos

% de 0 a 10
for n = 1 : Nt/10 -1

    vect_b = (M - 0.5*k*dt*R)*uhn + dt*M*(fi) - A*gi;
    vect_b01 = vect_b;
    vect_b01(fron_t) = 0;          

    % Resolver el sistema lineal
    whn = A01 \ vect_b01;
    uhn = whn + gi;
    
trisurf(elem,xi,yi,uhn)
zlim([0 50])
    title(n*dt)
     pause(0.1)
end

for n = Nt/10 : Nt

    vect_b = (M - 0.5*k*dt*R)*uhn + dt*M*(fi) - A*gi;
    vect_b02 = vect_b;
    vect_b02(fron1) = 0;          

    % Resolver el sistema lineal
    whn = A02 \ vect_b02;
    uhn = whn + gi;
    
trisurf(elem,xi,yi,uhn)
    title(n*dt)
zlim([0 50])
     pause(0.1)
end

%% EJERCICIO 10: convección + condiciones o homo

%% 10.a) crank nicolson

clear;

% Definir la función 
f = @(x,y,t) 0.*x + 0.*y +0.*t;
a = 0; 
k = 0.01;
b1 = 2;
b2 = 1;

    % MALLADO
load("malla_ej10.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = unique([e(1,:) e(2,:)]); %no repetidos
fron1 = find((xi+yi < 0.5) & (xi==0 | yi==0) );

% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.2;           
Nt = (tf - t0) / dt;    

% Matrices de masas y rigidez 
[R M] = assema(p,t,1,1,0);


% MATRICES DE CONVECCIÓN
Cx = sparse(length(xi), length(xi));
Cy=Cx;
% para hallar las dphis utilizamos pdecgrad que calcula
% el gradiente de una función escalar sobre un dominio 2D (o 3D) generado
% por pdetool

for i=1:length(xi) % columnas de C
    phi_i=0*xi'; 
    phi_i(i)=1; % genera la funcion base que vale uno en su nodo, 0 en el resto en cada iteración

    [phi_i_x, phi_i_y]=pdecgrad(p,t,'1',phi_i); % esa funcion devuelve dos vectores, con las derivadas

    [R, M, b]=assema(p,t,1,1,phi_i_x);
    Cx(:,i)=b;
    [R, M, b]=assema(p,t,1,1,phi_i_y);
    Cy(:,i)=b;
end

% FIN MATRICES CONV
A=(1+a*dt*0.5)*M+dt*k*0.5*R+b1*dt*0.5*Cx+b2*dt*0.5*Cy;

A0 = A;
A0(fron_d, :) = 0;     
A0(:, fron_d) = 0;    
for i = fron_d
    A0(i,i) = 1;
end
      
% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(fron1) = 1;               

% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 0*yi' + 0;
whn = uhn;

trisurf(elem,xi, yi, uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt
    fi0= f(xi',yi', (n-1)*dt);
    fi1 = f(xi',yi',n*dt);
    vect_b = ((1-a*dt*0.5)*M-k*dt*0.5*R-b1*dt*0.5*Cx-b2*dt*0.5*Cy)*uhn + dt*0.5*M*(fi0+fi1)-A*gi;
    vect_b0 = vect_b;
    vect_b0(fron_d) = 0;          

    % Resolver el sistema lineal
    whn = A0 \ vect_b0;
    uhn = whn + gi;
    
    trisurf(elem,xi,yi,uhn)
    title(n*dt)
pause
    
end

%% EJERCICIO 11 : convección no constante

clear;

% Definir la función 
f = @(x,y,t) 0.*x + 0.*y +0.*t;
a = 0; 
k = 0.01;
b1=@(x,y) -y; b2=@(x,y) x;


    % MALLADO
load("malla_ej11.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = unique([e(1,:) e(2,:)]); %no repetidos
fron1 = find((xi+yi < 0.5) & (xi==0 | yi==0) );

% Intervalo tiempo
t0 = 0; 
tf = 100;
dt = 0.2;           
Nt = (tf - t0) / dt;    

% Matrices de masas y rigidez 
[R M] = assema(p,t,1,1,0);


% MATRICES DE CONVECCIÓN
Cx = sparse(length(xi), length(xi));
Cy=Cx;
% para hallar las dphis utilizamos pdecgrad que calcula
% el gradiente de una función escalar sobre un dominio 2D (o 3D) generado
% por pdetool

for i=1:length(xi) %Columnas de C
    b1_i = b1(xi(i),yi(i));
    b2_i = b2(xi(i),yi(i));
    phi_i=0*xi'; 
    phi_i(i)=1; %Genera la funcion base que vale uno en su nodo, 0 en el resto
    [phi_i_x, phi_i_y]=pdecgrad(p,t,'1',phi_i); %Devuelve el gradiente de una función escalar sobre un dominio generado por pdetool
    phi_i_x = b1_i.*phi_i_x;
    phi_i_y = b2_i.*phi_i_y;
    [R, M, b]=assema(p,t,1,1,phi_i_x);
    Cx(:,i)=b;
    [R, M, b]=assema(p,t,1,1,phi_i_y);
    Cy(:,i)=b;
end
% FIN MATRICES CONV
A=(1+a*dt*0.5)*M+dt*k*0.5*R+dt*0.5*Cx+dt*0.5*Cy;

A0 = A;
A0(fron_d, :) = 0;     
A0(:, fron_d) = 0;    
for i = fron_d
    A0(i,i) = 1;
end
      
% función g para condiciones de contorno homogéneas
gi = 0 * xi'; 
gi(fron1) = 1;               

% Condición inicial y funcion evaluada en los nodos
uhn = exp((yi'.^2+(xi'-0.5).^2)/-0.01);
whn = uhn;

trisurf(elem,xi, yi, uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt
    fi0= f(xi',yi', (n-1)*dt);
    fi1 = f(xi',yi',n*dt);
    vect_b = ((1-a*dt*0.5)*M-k*dt*0.5*R-dt*0.5*Cx-dt*0.5*Cy)*uhn + dt*0.5*M*(fi0+fi1);
    vect_b0 = vect_b;
    vect_b0(fron_d) = 0;          

    % Resolver el sistema lineal
    whn = A0 \ vect_b0;
    uhn = whn + gi;
    
    trisurf(elem,xi,yi,uhn)
    title(n*dt)
        axis([-1 1 -1 1 0 0.4]);

pause(0.05)
    
end