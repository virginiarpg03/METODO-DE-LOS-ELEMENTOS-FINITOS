%% EJERCICIO 4

clear;

% Definir la función 
f = @(x,y,t) 0.*x + 0.*y +0.*t;
a = 0; 
k = 0.01;
b1 = 1;
b2 = 0;

    % MALLADO
    %lo he mallado muyy poco porque como hay matrices de convección mi
    %ordenador va mal

load("malla_ej4pec2.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

% ponemos las fronteras, la frontera dirichlet no es toda, solo es la
% frontera1 + frontera3
fron_t = unique([e(1,:) e(2,:)]); %no repetidos
fron2 = find(xi==5);
fron1 = find(xi==0);
fron3 = setdiff(fron_t, [fron1; fron2]);
fron_d= union(fron1,fron3) ;


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
B= (M - 0.5*dt*R - b1*dt*0.5*Cx - b2*dt*0.5*Cy)

A0 = A;
A0(fron_d, :) = 0;     
A0(:, fron_d) = 0;    
for i = fron_d
    A0(i,i) = 1;
end
      
% función g para condiciones de contorno homogéneas, ahora depende del
% tiempo por lo que la tengo que meter en el bucle
gi1 = 0*xi';
gi0 =gi1;
% Condición inicial y funcion evaluada en los nodos
uhn = 0 * xi' + 0*yi' + 0;
whn = uhn;

trisurf(elem,xi, yi, uhn);
pause;

% resolver problemas elípticos
for n = 1:Nt

    gi0(fron1) = 1- cos(2*pi*(n-1)*dt);
    gi1(fron1) = 1- cos(2*pi*n*dt);

    vect_b= gi0*(M - 0.5*dt*R - b1*dt*0.5*Cx - b2*dt*0.5*Cy) - gi1*A + whn*(M - 0.5*dt*R - b1*dt*0.5*Cx - b2*dt*0.5*Cy)
   % vect_b = gi0*B - gi1*A + whn*B; 
    vect_b0 = vect_b;
    vect_b0(fron_d) = 0;          

    % Resolver el sistema lineal
    whn = A0 \ vect_b0;
    uhn = whn + gi1;
    
    trisurf(elem,xi,yi,uhn)
    title(n*dt)
pause
    
end

%no entiendo porque da error al hacer el vector b si todas las variables
%tienen los mismos elementos

%% EJERCICIO 5

f1=@(x,y) 0*x;
f2=@(x,y) 0*x;

if 1

    myPDE=createpde;   %creación estructura de datos (para acceder myPDE. )
    
    geometryFromEdges(myPDE,g);   %creación mallado     g viene de pdetool y esport boundary
    %sirve para mallados 3d
    
    generateMesh(myPDE,'hmax',0.1,'geometricorder','quadratic'); %genera mallado (en 2d de triángulos) y lo guarda en estructura de datos
    
    
    xi=myPDE.Mesh.Nodes(1,:);
    yi=myPDE.Mesh.Nodes(2,:);
    elem=myPDE.Mesh.Elements(1:3,:)';

    applyBoundaryCondition(myPDE,'neumann','Edge',1:size(g,2),'q',1);
    
    specifyCoefficients(myPDE,'a',1,'c',1,'m',0,'d',0,'f',0);
    FEM_KM=assembleFEMatrices(myPDE,'none');   %estructura de datos
    M=FEM_KM.A;    %matriz de masas
    R=FEM_KM.K;    %matriz de rigidez
    
    [Cxl Cyl]=calcular_matrices_stokes(myPDE);    %matrices Cx y Cy
end

fron=find(diag(FEM_KM.Q)>0);   %nos devuelve los nodos de la frontera

fron_n=find(xi==5);
fron_d=setdiff(fron,fron_n)';
fron_d=[fron_d fron_d+length(xi)];
ind_1= find(xi==0)
ind_2=find(yi==2);

A=[R 0*R -Cxl; 0*R R -Cyl; Cxl' Cyl' zeros(size(Cxl,2))];

A0=A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end

f1i=f1(xi)';
f2i=f2(xi)';

g1i=0*xi';
g2i=0*xi';
g1i(ind_1)= yi(ind_1).*(0.5-yi(ind_1));
g2i(ind_2)= -9*(11/4 - xi(ind_2).*(xi(ind_2)-9/4));
gg=[g1i;g2i;zeros(size(Cxl,2),1)];

vect_b=[M*f1i; M*f2i; zeros(size(Cxl,2),1)];
vect_b=vect_b-A*gg;
vect_b(fron_d)=0;

wh=A0\vect_b;
uh=wh+gg;

u1h=uh(1:length(xi));
u2h=uh([1:length(xi)]+length(xi));
ph=uh((2*length(xi)+1):end);

figure
trisurf(elem,xi,yi,u1h)

figure
trisurf(elem,xi(1:size(Cxl,2)),yi(1:size(Cyl,2)),ph)    %dibujo de la presión

figure
quiver(xi,yi,u1h',u2h',4),axis equal