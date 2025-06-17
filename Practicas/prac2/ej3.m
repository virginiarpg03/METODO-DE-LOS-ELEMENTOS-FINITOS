%% ejercicio 3

% hacemos mallado con pdetool
load("malla_ej3.mat")
f = @(x,y) 0*x + 1;

a=0;
k=1;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';
fron_d = unique([e(1,:) e(2,:)]); %no repetidos

%buscar fronteras
ind_4 = find(yi==0);
ind_3 = find(xi==0);
ind_2 = find(abs(xi.^2 + yi.^2 -1 < 0.001));  % a veces si lo ponemos =0 da error
ind_1 = setdiff(fron_d,[ind_4 ind_2 ind_3]); %diferencia de conjuntos


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
vect_b = M*fi ; 
vect_b(fron_d) = 0;

uh = A0\vect_b;


trisurf(elem,xi,yi,uh)



%%
%b)

f = @(x,y) 0*x;

fi = f(xi,yi)';



gi = 0*xi' ; %puede valer cualquier cosa dentro, pues que valga cero
 %mas condiciones para la frontera

gi(ind_1)=0;
gi(ind_2)=1;
gi(ind_3)=tanh(2*(yi(ind_3)+3));

vect_b = M*fi - A*gi; 
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh2 = wh + gi;

figure 

trisurf(elem,xi,yi,uh2)


%% 