%% ejercicio 7


f = @(x,y) 1 + 0.1*(x + y) ;
load("malla_ej7.mat")
a=1;
k=1;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';
fron_d = unique([e(1,:) e(2,:)]); %no repetidos

%buscar fronteras a cachos
ind_1 = find(abs(xi.^2 + yi.^2 -4 == 0));
ind_2 = find(yi==0);
ind_3 = find(xi==0);

ind_4 = setdiff(fron_d,[ind_1 ind_2 ind_3]); %diferencia de conjuntos


%%
% comprobar fronteras para ver si lo he hecho bien
figure
hold on
scatter(xi(fron_d),yi(fron_d),'b')
scatter(xi(ind_1),yi(ind_1),'r','.')
hold off
%%
[R M]=assema(p,t,1,1,0);

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
gi(ind_2)=2;
gi(ind_3)=2-xi;


vect_b = M*fi - A*gi; % caract del problema, viene a ser L() paratod j
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh = wh + gi;

trisurf(elem,xi,yi,uh)