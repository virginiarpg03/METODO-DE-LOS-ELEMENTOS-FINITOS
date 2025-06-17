%% mallado de pdetool y cuadratura con los vértices


load('malla_ej9.mat')

f=@(x,y) exp(-((x-1).^2 +y.^2)/0.4^2);

xi=p(1,:);    %valor de la x de cada nodo
yi=p(2,:);    %valor de la y de cada nodo

elem=t(1:3,:)';    %guarda los nodos de cada triángulo


trisurf(elem,xi,yi,f(xi,yi)) %dibujo


Ne=length(elem);    %nos da elnúmero de triángulos que ha generado el mallado

Xg=[0 1 0;0 0 1];    %los valores de (xg,yg) que se usan para la cuadratura
wi=[1/6 1/6 1/6];    %pesos que se usan en la cuadratura
Area=zeros(Ne,1);

Q=0;

for i=1:Ne    %las coordenadas de cada nodo de cada triángulo
    X1=p(:,elem(i,1));
    X2=p(:,elem(i,2));
    X3=p(:,elem(i,3));

%ahora aplicamos la transformación afín para triángulos 
bi=X1;
Ai=[(X2-X1) (X3-X1)];
Area(i)=0.5*abs(det(Ai));

    for j=1:3
        X=bi+Ai*Xg(:,j); %creamos los tres vértices transformados 
        Q=Q+wi(j)*abs(det(Ai))*f(X(1), X(2));
     end 
end 
Q
A=sum(Area)
