f1=@(x,y) 0*x;
f2=@(x,y) 0*x;

if 0

    myPDE=createpde;   %creación estructura de datos (para acceder myPDE. )
    
    geometryFromEdges(myPDE,g);   %creación mallado     g viene de pdetool y esport boundary
    %sirve para mallados 3d
    
    generateMesh(myPDE,'hmax',0.1,'geometricorder','quadratic'); %genera mallado (en 2d de triángulos) y lo guarda en estructura de datos
    
    
    xi=myPDE.Mesh.Nodes(1,:);
    yi=myPDE.Mesh.Nodes(2,:);
    elem=myPDE.Mesh.Elements(1:3,:)';
    
    specifyCoefficients(myPDE,'a',1,'c',1,'m',0,'d',0,'f',0);
    FEM_KM=assembleFEMatrices(myPDE,'none');   %estructura de datos
    M=FEM_KM.A;    %matriz de masas
    R=FEM_KM.K;    %matriz de rigidez
    
    [Cxl Cyl]=calcular_matrices_stokes(myPDE);    %matrices Cx y Cy
end

fron_d=find(diag(FEM_KM.Q)>0);
fron_d=[fron_d fron_d+length(xi)];
ind_1=find(xi==3.25 & yi>=0 & yi<=0.08 | xi==4.75 & yi>=0 & yi<=0.08 | yi==0.08 & xi>=3.25 & xi<=4.75);
ind_2=find(yi==3 & xi>=0 & xi<=0.1 | yi==4 & xi>=0 & xi<=0.1 | xi==0.1 & yi>=3 & yi<=4);

A=[R 0*R -Cxl; 0*R R -Cyl; Cxl' Cyl' zeros(size(Cxl,2))];

A0=A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end

f1i=f1(xi)';
f2i=f2(xi)';

g1i=0*xi';
g2i=0*xi';
g1i(ind_2)=3;   %hay que elegir una
gg=[g1i;g2i;zeros(size(Cxl,2),1)];

vect_b=[M*f1i; M*f2i; zeros(size(Cxl,2),1)];
vect_b=vect_b-A*gg;
vect_b(fron_d)=0;

wh=A0\vect_b;
bh=wh+gg;

b1h=bh(1:length(xi));
b2h=bh([1:length(xi)]+length(xi));
ph=bh((2*length(xi)+1):end);

%trisurf(elem,xi,yi,b1h)

trisurf(elem,xi(1:size(Cxl,2)),yi(1:size(Cyl,2)),ph)    %dibujo de la presión

%quiver(xi,yi,u1h',u2h',4),axis equal