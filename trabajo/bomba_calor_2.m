%Ecuación del calor

%Trabajo, parte bomba de calor

bomba_calor

f=@(x,y,t) 0*x;

a=0;
k=0.025; 

t0=0;
tf=100;
dt=0.1;
Nt=(tf-t0)/dt;

xi=p(1,:);
yi=p(2,:);
elem=t(1:3,:)';

front_d=unique([e(1,:) e(2,:)]);
ind_1=find(xi==3.25 & yi>=0 & yi<=0.08 | xi==4.75 & yi>=0 & yi<=0.08 | yi==0.08 & xi>=3.25 & xi<=4.75);
ind_2=find(yi==3 & xi>=0 & xi<=0.1 | yi==4 & xi>=0 & xi<=0.1 | xi==0.1 & yi>=3 & yi<=4);
ind_3=find(xi==8 & yi>=6.85 & yi<=7.15);


% Crear matrices de convección
Cx_b1 = sparse(length(xi), length(xi));
Cy_b2 = sparse(length(xi), length(xi));

for i = 1:length(xi)

    % Inicializamos phi_i para el nodo i-ésimo
    phi_i = 0 * xi'; 
    phi_i(i) = 1; % phi_i es 1 en el nodo i y 0 en el resto
    
    % Gradientes de la función phi_i
    [phi_i_x, phi_i_y] = pdecgrad(p, t, 1, phi_i);
    
    % Incorporar b1h (velocidad en x) en la matriz de convección
    b1_local = b1h(i); % Valor de b1h en el nodo i
    [R, M, b] = assema(p, t, 1, 1, phi_i_x * b1_local); % Multiplicamos gradiente por b1h
    Cx_b1(:, i) = b;

    % Incorporar b2h (velocidad en y) en la matriz de convección
    b2_local = b2h(i); % Valor de b2h en el nodo i
    [R, M, b] = assema(p, t, 1, 1, phi_i_y * b2_local); % Multiplicamos gradiente por b2h
    Cy_b2(:, i) = b;

end


%MATRIZ A
A = M + dt*0.5*Cx_b1 + dt*0.5*Cy_b2 + k*dt*0.5*R;
A0=A;
A0(front_d,:)=0;
A0(:,front_d)=0;
for i=front_d
    A0(i,i)=1;
end


gi=0*xi';
gi(ind_1)=5;
gi(ind_2)=30; 
gi(ind_3)=17;

uhn=0*xi'+20;

trisurf(elem,xi,yi,uhn)
title(0)
pause
for n=1:Nt 
   vect_b = (M - dt*0.5*Cx_b1 - dt*0.5*Cy_b2 - k*dt*0.5*R)*uhn - A*gi;
    vect_b(front_d) = 0;

    whn = A0\vect_b;
    uhn = whn + gi;

    trisurf(elem,xi,yi,uhn)
    title(n*dt)
    view(2)
    shading interp
    colorbar
    pause
end
%En el primer paso de tiempo tengo que observar las condiciones de
%contorno.


