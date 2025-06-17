% --- Implementación resolución numérica PEF elíptico en 2D segundo ejemplo---

load("malla_ej9.mat")

f = @(x,y) 0*x + 2;

a=0;
k=1;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';

%esto falla muy a menudo por la comparación, ponerle una tolerancia, 
fron_d = find(abs(xi.^2/4+yi.^2/16-1)<0.001);

[R M]=assema(p,t,1,1,0);

A = a*M + k*R;
%¿como eliminar para sacar A0?
A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    



fi = f(xi,yi)';

vect_b = M*fi; %en general asi, 
vect_b(fron_d) = 0;

uh = A0\vect_b;

trisurf(elem,xi,yi,uh)
%para analizar, pintar puntos de frontera dirichlet 
%plot(xi(fron_d),)