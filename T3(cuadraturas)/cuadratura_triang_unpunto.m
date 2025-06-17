%% utilizando mallado de pdetool y un solo nodo de integración

f=@(x,y) 1+exp(-0.5*(x.^2 + y.^2)); 

elem=t(1:3,:)';
xi=p(1,:);
yi=p(2,:);

Ne=length(t);

Xg=[1/3; 1/3];

Qf=0;
area_total=0;

for i=1:Ne
    X1=p(:,t(1,i));  %primer nodo elemento i
    X2=p(:,t(2,i));
    X3=p(:,t(3,i));

    A= [X2-X1 X3-X1];

    X=X1+A*Xg;

    area=0.5*det(A);
    Qf= Qf+area*f(X(1),X(2)); %cuadratura en el elemento i-ésimo
    area_total=area_total+area;
end

Qf
area_total
%para dibujar la superficie: trisurf(elem,xi,yi,f(xi,yi))
%si voy refinando el mallado (haciendo más triángulos) debe ir mejorando la
%aproximación