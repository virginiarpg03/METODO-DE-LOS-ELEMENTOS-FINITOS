
xi=p(1,:);
yi=p(2,:);

%%ECUACIÓN MATRICIAL
[R,M]=assema(p,t,1,1,0);

Cx = sparse(length(xi), length(xi));Cy=Cx;
for i=1:length(xi) %Columnas de C
    phi_i=0*xi'; 
    phi_i(i)=1; %Genera la funcion base que vale uno en su nodo, 0 en el resto
    [phi_i_x, phi_i_y]=pdecgrad(p,t,'1',phi_i); %Devuelve el gradiente de una función escalar sobre un dominio generado por pdetool
    [R, M, b]=assema(p,t,1,1,phi_i_x);
    Cx(:,i)=b;
    [R, M, b]=assema(p,t,1,1,phi_i_y);
    Cy(:,i)=b;
end

A=a*M+k*R+b1*Cx+b2*Cy;

fi=f(xi,yi);
gi=0*xi;
frontera=unique([e(1,:) e(2,:)]);
ind_1=find(xi+yi<0.5);
ind_2=find(xi+yi>0.5);
gi(fron1)=1;
gi(ind_2)=0;
vect_b=M*fi'-A*gi';

%Ahora aplico las condiciones de contorno para construir A0 y vect_b0
frontera_d=frontera;                     
A0=A;
A0(:,frontera_d)=0; 
A0(frontera_d,:)=0; 

for i=frontera_d
    A0(i,i)=1; 
end

vect_b0=vect_b;
vect_b0(frontera_d)=0;

%SOLUCIÓN Y GRÁFICO
wh=A0\vect_b0; 
uh=wh+gi';
figure
trisurf(t(1:3,:)',xi,yi,uh)