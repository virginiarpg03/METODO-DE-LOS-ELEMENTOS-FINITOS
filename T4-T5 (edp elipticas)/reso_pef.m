% --- Implementación resolución numérica PEF elíptico en 1D ---

% Sea PEF con condiciones dirichlet homogéneas, u-u''=f , solución exacta 
% u = sen(nx), f = (1+pi^2)*sen(nx)

f = @(x) 0*x+1;

a=1;
k=0.1;

c=0;
d=1;

Ne = 20;
h = (d-c)/Ne;
xi = c:h:d;


matriz_masa_rigidez_elemref_lineal
A = a*M + k*R;

A0 = A;
A0(1,:)=0; %se puede compactar A0([1 end],:)=0;
A0(:,1)=0;
A0(end,:)=0;
A0(:,end)=0;

A0(1,1)=1;
A0(end,end)=1;

fi = f(xi');
vect_b = M*fi ; 
vect_b([1 end]) = 0;

uh = A0\vect_b;
plot(xi,uh)
hold on
%plot(xi_pint,sol(xi_pint), 'g')
hold off