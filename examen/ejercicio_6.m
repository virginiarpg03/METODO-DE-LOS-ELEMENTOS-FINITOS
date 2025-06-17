%% ejercicio 6

% hacemos el malado con pdetool
load("malla_ej6_examen.mat")

f=@(x,y) (x.*exp(x-y.^2/9) .*sin(y))    ./     (x +y +1); 

elem=t(1:3,:)';
xi=p(1,:);
yi=p(2,:);

Ne=length(t);

Xg=[1/3; 1/3];

Qf=0;

for i=1:Ne
    X1=p(:,t(1,i));  %primer nodo elemento i
    X2=p(:,t(2,i));
    X3=p(:,t(3,i));

    A= [X2-X1 X3-X1];

    X=X1+A*Xg;

    area=0.5*det(A);
    Qf= Qf+area*f(X(1),X(2)); %cuadratura en el elemento i-ésimo
end

Qf
%dubijo
trisurf(elem,xi,yi,f(xi,yi))


%% apartado b)

% lo hago con la cuadratura de grado 1 w=1, x= (1+- raiz(2)) /12
f=@(x,y) (x.*exp(x-y.^2/9) .*sin(y))    ./     (x +y +1); 

load("malla_ej6_examen.mat")

xi=p(1,:);    %valor de la x de cada nodo
yi=p(2,:);    %valor de la y de cada nodo

elem=t(1:3,:)';    %guarda los nodos de cada triángulo

Ne=length(elem);    %nos da elnúmero de triángulos que ha generado el mallado

Xg=[ (1-sqrt(2))/12    ; (1-sqrt(2))/12 ];    
wi=[1];    

Q=0;

for i=1:Ne    %las coordenadas de cada nodo de cada triángulo
    X1=p(:,elem(i,1));
    X2=p(:,elem(i,2));
    X3=p(:,elem(i,3));

%ahora aplicamos la transformación afín para triángulos 
bi=X1;
Ai=[(X2-X1) (X3-X1)];

   
        X=bi+Ai*Xg(:,1); %creamos los tres vértices transformados 
        Q=Q+wi(1)*abs(det(Ai))*f(X(1), X(2));
      
end 
Q
trisurf(elem,xi,yi,f(xi,yi))


%%

