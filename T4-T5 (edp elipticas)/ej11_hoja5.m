% --- Implementación resolución numérica PEF elíptico en 2D segundo ejemplo---

%ej11 (15) 

f = @(x,y) 0*x;

a=0;
k=1;

b1 = 2;
b2 = 1;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';
fron_d = find(xi.*yi==0 | xi == 1 | yi == 1); 
%fron_d = unique(e(1,:)) %una de dos
ind_1 = find((xi==0 & yi<0.5) |(yi==0 & xi<0.5));


calcular_matrices_convec_masa_rig
A = b1*Cx + b2*Cy+k*R;
%¿como eliminar para sacar A0?
A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    



%fi = f(xi,yi)';
gi = 0*xi';
gi(ind_1)=1;


%espabilar y darse cuenta de que estas definiendo valores en los nodos, g
%claro que esta dentro del esp de elemfin y es contínua


vect_b = -A*gi; %en general asi, 
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh = wh + gi;

trisurf(elem,xi,yi,uh)
%shading interp