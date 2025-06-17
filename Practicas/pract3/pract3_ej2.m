%Práctica 3, Ejercicio 2 apartado A

f1=@(x,y) 0*x;
f2=@(x,y) 0*x;

if 1

    myPDE=createpde;   %creación estructura de datos (para acceder myPDE. )
    
    geometryFromEdges(myPDE,g);   %creación mallado     g viene de pdetool y esport boundary
    %sirve para mallados 3d
    
    generateMesh(myPDE,'hmax',0.1,'geometricorder','quadratic'); %genera mallado (en 2d de triángulos) y lo guarda en estructura de datos
    
    
    xi=myPDE.Mesh.Nodes(1,:);
    yi=myPDE.Mesh.Nodes(2,:);
    elem=myPDE.Mesh.Elements(1:3,:)';

    applyBoundaryCondition(myPDE,'neumann','Edge',1:size(g,2),'q',1);
    
    specifyCoefficients(myPDE,'a',1,'c',1,'m',0,'d',0,'f',0);
    FEM_KM=assembleFEMatrices(myPDE,'none');   %estructura de datos
    M=FEM_KM.A;    %matriz de masas
    R=FEM_KM.K;    %matriz de rigidez
    
    [Cxl Cyl]=calcular_matrices_stokes(myPDE);    %matrices Cx y Cy
end

fron=find(diag(FEM_KM.Q)>0);   %nos devuelve los nodos de la frontera
fron_n=find(yi==5 & xi>=8.5 & yi<=9 | xi==15 & yi>=-0.2 & yi<=0);
fron_d=setdiff(fron,fron_n)';
fron_d=[fron_d fron_d+length(xi)];

A=[R 0*R -Cxl; 0*R R -Cyl; Cxl' Cyl' zeros(size(Cxl,2))];

A0=A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end

f1i=f1(xi)';
f2i=f2(xi)';

ind_2=find(xi==0 & yi>=-1 & yi<=1);
g1i=0*xi';
g2i=0*xi';
g1i(ind_2)=1-yi(ind_2).^2;
gg=[g1i;g2i;zeros(size(Cxl,2),1)];

vect_b=[M*f1i; M*f2i; zeros(size(Cxl,2),1)];
vect_b=vect_b-A*gg;
vect_b(fron_d)=0;

wh=A0\vect_b;
uh=wh+gg;

u1h=uh(1:length(xi));
u2h=uh([1:length(xi)]+length(xi));
ph=uh((2*length(xi)+1):end);

figure
trisurf(elem,xi,yi,u1h)

figure
trisurf(elem,xi(1:size(Cxl,2)),yi(1:size(Cyl,2)),ph)    %dibujo de la presión

figure
quiver(xi,yi,u1h',u2h',4),axis equal