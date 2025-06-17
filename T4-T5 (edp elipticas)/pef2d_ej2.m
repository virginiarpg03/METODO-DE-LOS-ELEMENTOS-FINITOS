% --- Implementación resolución numérica PEF elíptico en 2D segundo ejemplo---


f = @(x,y) 0*x;

a=1;
k=1e-2;

xi = p(1,:); %del mallado importado
yi = p(2,:);
elem = t(1:3,:)';
fron_d = find(yi==1 | xi == 1); %no unique pq solo queremos un trozo cabu

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
gi = 0*xi' ; %puede valer cualquier cosa dentro, pues que valga cero
%mas condiciones para la frontera
gi(fron_d)=4*xi(fron_d) + 4;
%para ind_2 ya era cero


vect_b = M*fi - A*gi; %en general asi, 
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh = wh + gi;

trisurf(elem,xi,yi,uh)