load("malladoej52022.mat", "p","e", "t") 

f=@(x,y) 1 + exp(-0.5*(x.^2 + y.^2));
xi=p(1,:);
yi=p(2,:);
elem=t(1:3,:)'; %matriz que tiene en las tres primeras filas los nodos que forman el triángulo
Ne=length(elem); %Ne es el número de triángulos que tengo

w1=1/2; %peso para cuadratura de gauss legendre en triángulos de n=1
V1g=[1/3;1/3]; 

Q22=0;

for i=1:Ne
    %Sacamos los nodos del elemento i
    X1=p(:,elem(i,1));
    X2=p(:,elem(i,2));
    X3=p(:,elem(i,3));
    
    %Transfromación afín X=bi+Ai*X1g
    bi=X1;
    Ai=[(X2-X1) (X3-X1)]; 

    V1=bi+Ai*V1g;
    Q22= Q22+ abs(det(Ai))*(w1*f(V1(1),V1(2)));
end 

masa= Q22*2e-3