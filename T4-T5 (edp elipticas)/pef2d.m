% --- Implementación resolución numérica PEF elíptico en 2D ---


f = @(x,y) 4000*exp(-10*(x.^2+y.^2));

a=0;
k=1;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';
fron_d = unique([e(1,:) e(2,:)]); %no repetidos

%buscar fronteras a cachos
ind_1 = find(yi==0);
ind_2 = find(xi==0);
ind_3 = setdiff(fron_d,[ind_1 ind_2]); %diferencia de conjuntos




[R M]=assema(p,t,1,1,0);
% los 1,1,0 son a, k, 0 que matlab las pone dentro, mejor fuera, el cero ya veremos
%proceso laborioso, sencillo

A = a*M + k*R; % caract del problema, viene a ser A() paratod j
%¿como eliminar para sacar A0?
A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    



fi = f(xi,yi)';
gi = 0*xi' ; %puede valer cualquier cosa dentro, pues que valga cero
%mas condiciones para la frontera
gi(ind_1)=12*xi(ind_1);
gi(ind_3)=4*xi(ind_3).^2;
%para ind_2 ya era cero


vect_b = M*fi - A*gi; % caract del problema, viene a ser L() paratod j
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh = wh + gi;

trisurf(elem,xi,yi,uh)