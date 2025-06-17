%para implementar funciones base cuadráticas, usamos la estructura de problema para PDE

f1 = @(x,y) 0*x;
f2 = @(x,y) 0*x;

if 1 %para que esta parte no se ejecute (como es muy lenta solo la ejecuto 1 vez)
    myPDE = createpde;
    
    %con pdetool le doy a export geometry y me saca g,b
    geometryFromEdges(myPDE, g);
    %también admite geometrías es .stl, por tanto esto nos vale para un problema en 3D
    
    generateMesh(myPDE,'hmax',0.1,'geometricorder','quadratic');
    
    xi = myPDE.Mesh.Nodes(1,:);
    yi = myPDE.Mesh.Nodes(2,:);
    elem = myPDE.Mesh.Elements(1:3,:)'; 
    %por ser cuadraticos, tenemos 6 nodos por triangulo (vértices y puntos medios de arista)
    %solo cogemos los vertices porque los usaremos para el plot
    
    applyBoundaryCondition(myPDE,'neumann','Edge',1:size(g,2),'q',1);
    specifyCoefficients(myPDE,'a',1,'c',1,'m',0,'d',0,'f',0);
    FEM_KM = assembleFEMatrices(myPDE,'none'); %obtenemos una estructura de datos
    M = FEM_KM.A;
    R = FEM_KM.K;
    
    [Cxl Cyl] = calcular_matrices_stokes(myPDE);
end

fron = find(diag(FEM_KM.Q)>0); %solo los nodos de la frontera
fron_n = find(xi==5 & abs(yi)<0.5);
fron_d = setdiff(fron,fron_n)';
fron_d = [fron_d fron_d+length(xi)]; %también hay que anular los nodos para u2

A = [R 0*R -Cxl; 0*R R -Cyl; Cxl' Cyl' zeros(size(Cxl,2))];

A0 = A;
A0(fron_d,:) = 0;
A0(:,fron_d) = 0;

for i=fron_d
    A0(i,i) = 1;
end

f1i = f1(xi)';
f2i = f2(xi)';

ind = find(xi==0);
g1i = 0*xi';
g2i = 0*xi';
g1i(ind) = 0.5^2 - yi(ind).^2;
gg = [g1i; g2i; zeros(size(Cxl,2),1)];

vect_b = [M*f1i; M*f2i; zeros(size(Cxl,2),1)];
vect_b = vect_b - A*gg;
vect_b(fron_d) = 0;

wh = A0\vect_b;
uh = wh + gg;

u1h = uh(1:length(xi));
u2h = uh([1:length(xi)]+length(xi));
ph = uh((2*length(xi)+1):end);

figure(1)
trisurf(elem,xi,yi,u1h)

figure(2)
trisurf(elem,xi,yi,u2h) %vemos error de máquina. debería ser 0

figure(3)
quiver(xi,yi,u1h',u2h',4), axis equal

figure(4)
trisurf(elem,xi(1:size(Cxl,2)),yi(1:size(Cxl,2)),ph) %solo numero de nodos lineales