function [M,R]=MRlineales(xi,Ne)
    %Funciones base y sus derivadas
    phi0_g = @(x) 0.5*(1-x); 
    phi1_g = @(x) 0.5*(1+x);
    dphi0_g = @(x) 0*x -0.5; 
    dphi1_g = @(x) 0*x + 0.5;
    
    %Parámetros de cuadratura Gauss Legrende para 2 puntos
    gi=[1 1]; 
    xgi=[-1/sqrt(3) 1/sqrt(3)/3]; 
    
    %Relleno matriz de masas y rigidez de referencia:Usamos sparce porque son matrices "vacías" en vez de llenas de ceros y son más eficientes
    M_g = sparse(2,2);
    M_g(1,1) = sum(gi.*phi0_g(xgi).*phi0_g(xgi));
    M_g(1,2) = sum(gi.*phi0_g(xgi).*phi1_g(xgi));
    M_g(2,1) = sum(gi.*phi1_g(xgi).*phi0_g(xgi));
    M_g(2,2) = sum(gi.*phi1_g(xgi).*phi1_g(xgi));
    R_g = sparse(2,2);
    R_g(1,1) = sum(gi.*dphi0_g(xgi).*dphi0_g(xgi));
    R_g(1,2) = sum(gi.*dphi0_g(xgi).*dphi1_g(xgi));
    R_g(2,1) = sum(gi.*dphi1_g(xgi).*dphi0_g(xgi));
    R_g(2,2) = sum(gi.*dphi1_g(xgi).*dphi1_g(xgi));
    
    %Relleno las matrices de masas y rigidez:Al hacer operaciones con matrices sparce solo se alamcenan los resultados no nulos
    Nt = length(xi);
    M = sparse(Nt,Nt);
    R = sparse(Nt,Nt);
    
    %En los elementos de la diagonal principal se van solapando los valores de las matrices de referencia 
    for i=1:Ne
        hi=xi(i+1)-xi(i); %Por si hubiera usado nodos no equiespaciados
        M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+hi*0.5*M_g;
        R(i:i+1,i:i+1) = R(i:i+1,i:i+1)+(2/hi)*R_g;
    end

end