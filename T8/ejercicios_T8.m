%% EJERCICIOS ELASTOESTÁTICA

%% EJERCICIO 6 (barra condiciones homo)

    % PARÁMETROS
E = 2.7e8;
nu = 0.2;
rho = 1e3;
g = 9.8;

mu = E*0.5/(1+nu);
lambda = E*nu/((1+nu)*(1-2*mu));


    % MALLADO
load("mallado_barra.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = find( (xi<=0.5 & yi==0) | (xi<=0.5 & yi ==0.2) | (xi>=10.5 & yi==0) | (xi>=10.5 & yi==0.2) | xi==0 | xi==11);
fron_d = [fron_d fron_d+length(xi)];    % para que afecte tanto a uh1 como uh2


    % CREAR MATRICES R
[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R12 M] = assema(p,t,[0 1 0 0]',1,0);
[R21 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);


    % DIVISIÓN DE LA MATRIZ A
A11 = (lambda+2*mu)*R11 + mu*R22;
A12 = lambda*R21 + mu*R12;
A21 = lambda*R12 + mu*R21;
A22 = (lambda+2*mu)*R22 + mu*R11;

A = [A11 A12; A21 A22];

    %VECTOR b
f1i = 0*xi';
f2i = 0*xi' - rho*g;

vect_b1 = M*f1i;
vect_b2 = M*f2i;

vect_b = [vect_b1; vect_b2];

    % CONDICIONES DE CONTORNO
A0 = A;
A0(fron_d,:) = 0;
A0(:,fron_d) = 0;
for i = fron_d
    A0(i,i) = 1;
end
vect_b(fron_d) = 0;

    % SOLUCIÓN
uh = A0\vect_b;
u1h = uh(1:length(xi));
u2h = uh(length(xi)+1:end);


figure(1)                       % deformacion direccion x
trisurf(elem,xi,yi,u1h)
axis equal
figure(2)                       % deformacion direccion y
trisurf(elem,xi,yi,u2h)
figure(3)                       % deformacion direccion y vista como la barra
trisurf(elem,xi,yi,u2h)
view(2)
shading interp
colorbar
axis equal
figure(4)                       % deformación total
trisurf(elem,xi+u1h',yi+u2h',0*xi)
axis equal


%% EJERCICIO 8
%save("malla_ej8","t","e","p")

    % PARÁMETROS
E = 2.7e8;
nu = 0.2;

mu = E*0.5/(1+nu);
lambda = E*nu/((1+nu)*(1-2*mu));

    % MALLADO
load("malla_ej8.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = find( yi==0 | yi==4); %La frontera dirichlet es donde tenemos desplazamiento fijado
ind1= find(yi==0);
ind2= find(yi==4);

fron_d=[fron_d fron_d+length(xi)];

  % CREAR MATRICES R
[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R12 M] = assema(p,t,[0 1 0 0]',1,0);
[R21 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);


 % DIVISIÓN DE LA MATRIZ A
A11 = (lambda+2*mu)*R11 + mu*R22;
A12 = lambda*R21 + mu*R12;
A21 = lambda*R12 + mu*R21;
A22 = (lambda+2*mu)*R22 + mu*R11;

A = [A11 A12; A21 A22];

A0 = A;
A0(fron_d,:) = 0;
A0(:,fron_d) = 0;
for i = fron_d
    A0(i,i) = 1;
end


    %VECTOR b
f1i = 0*xi';
f2i = 0*xi';
g1i = 0*xi';
g2i = 0*xi';
g1i(ind1) = 0;
g2i(ind2) = 0.5;
gi= [g1i ; g2i];

vect_b1 = M*f1i ;
vect_b2 = M*f2i ;
vect_b = [vect_b1; vect_b2];

vect_b = vect_b - A*gi;

vect_b(fron_d) = 0;


    % RESOLVER SISTEMA
wh=A0\vect_b;
uh=wh+gi;
u1h=uh(1:length(xi));
u2h=uh(length(xi)+1:end);


figure(1)
trisurf(elem,xi,yi,u1h) %Deformación x

figure(2)
trisurf(elem,xi,yi,u2h) %Deformación y

figure(3)
trisurf(elem,xi,yi,u2h) %Deformación y vista desde arriba con colorbar
view(2)
shading interp%Interesa para ver deformaciones pequeñas
colorbar
axis equal

figure(4)
trisurf(elem,xi+u1h',yi+u2h',0*xi) %Barra deformada
view(2)
axis equal


%% ejercicio random con newman no homo

   %pdetool -> Pinto dominio, le doy a triangular (3/4 veces) y exporto el mesh
%save('ej6.mat', 'p', 'e', 't') -> guardo el mallado en un .mat
load("malla_ej8.mat", "p","e", "t") 
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';
fron_d = find(yi==0); %La frontera dirichlet es donde tenemos desplazamiento fijado
ind2 = find(yi==4);

[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R12 M] = assema(p,t,[0 1 0 0]',1,0);
[R21 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);

E = 2.7e8;
nu = 0.2;
mu = E*0.5/(1+nu);
lambda = E*nu/((1+nu)*(1-2*mu));

A11 = (lambda+2*mu)*R11 + mu*R22;
A12 = lambda*R21 + mu*R12;
A21 = lambda*R12 + mu*R21;
A22 = (lambda+2*mu)*R22 + mu*R11;
A = [A11 A12; A21 A22];

f1i = 0*xi';
f2i = 0*yi';
gN1i= 0*xi';
gN2i= 0*yi';
gN2i(ind2)=-5;

vect_b1 = M*f1i + M*gN1i;
vect_b2 = M*f2i + M*gN2i;
vect_b = [vect_b1; vect_b2];

%Aplico condiciones de contorno
fron_dt= [fron_d fron_d+length(xi)]; 
A0 = A;
A0(fron_dt,:) = 0;
A0(:,fron_dt) = 0;
for i = fron_dt
    A0(i,i) = 1;
end
vect_b0=vect_b;
vect_b0(fron_dt) = 0;

%Resuelvo el sistema
uh = A0\vect_b0;
u1h = uh(1:length(xi));
u2h = uh(length(xi)+1:end);

figure(1)
trisurf(elem,xi,yi,u1h) %Deformación x

figure(2)
trisurf(elem,xi,yi,u2h) %Deformación y

figure(3)
trisurf(elem,xi,yi,u2h) %Deformación y vista desde arriba con colorbar
view(2)
shading interp%Interesa para ver deformaciones pequeñas
colorbar
axis equal

figure(4)
trisurf(elem,xi+u1h',yi+u2h',0*xi) %Barra deformada
view(2)
axis equal





