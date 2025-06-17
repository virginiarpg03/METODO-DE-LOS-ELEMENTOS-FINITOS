%% EJERCICIO 1    

% PARÁMETROS
mu = 10^8;
lambda = 7.5*10^7;
theta = pi/4;

    % MALLADO
load("malla_ejercicio1.mat")
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = find(yi==1 | yi==-1); %La frontera dirichlet es donde tenemos desplazamiento fijado
ind_1= find(yi==1);
ind_2= find(yi==-1);

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

g11i = 0*xi';
g12i = 0*xi';
g11i(ind_1) = (1-xi(ind_1)).*(1-cos(theta));
g12i(ind_1) = (1-xi(ind_1)).*sin(theta);
g1i= [g11i ; g12i];

g21i = 0*xi';
g22i = 0*xi';
g21i(ind_2) = (1-xi(ind_2)).*(1-cos(theta));
g22i(ind_2) = (xi(ind_2)-1).*sin(theta);
g2i= [g21i ; g22i];


vect_b1 = M*f1i ;
vect_b2 = M*f2i ;
vect_b = [vect_b1; vect_b2];

vect_b = vect_b - A*g1i -A*g2i;

vect_b(fron_d) = 0;


    % RESOLVER SISTEMA
wh=A0\vect_b;
uh=wh + g1i + g2i;
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
trisurf(elem,xi+u1h',yi+u2h',0*xi)    %Barra deformada
view(2)
axis equal