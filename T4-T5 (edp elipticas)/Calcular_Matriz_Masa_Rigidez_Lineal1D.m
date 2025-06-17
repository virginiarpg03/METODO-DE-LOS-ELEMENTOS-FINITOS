%% CREAR MATRICES DE MASA Y RIGIDEZ EN EL ELEMENTO DE REFERENCIA 1D (con cuadratura gauss legendre) Y ENSAMBLARLAS


%funciones base y sus derivadas
    phi0_g = @(x) 0.5*(1-x);
    phi1_g = @(x) 0.5*(1+x);
    
    dphi0_g = @(x) 0*x-0.5; 
    dphi1_g = @(x) 0*x + 0.5;
    
    % parámetros de cuadratura
    xgi = [-sqrt(3)/3 sqrt(3)/3]; % Nodos y pesos GL
    wi = [1 1];
    
    M_g = sparse(2,2);
    M_g(1,1) = sum(wi.*phi0_g(xgi).*phi0_g(xgi));
    M_g(1,2) = sum(wi.*phi0_g(xgi).*phi1_g(xgi));
    M_g(2,1) = sum(wi.*phi1_g(xgi).*phi0_g(xgi));
    M_g(2,2) = sum(wi.*phi1_g(xgi).*phi1_g(xgi));
    
    % MATRIZ DE RIGIDEZ DE REFERENCIA
    R_g = sparse(2,2);
    R_g(1,1) = sum(wi.*dphi0_g(xgi).*dphi0_g(xgi));
    R_g(1,2) = sum(wi.*dphi0_g(xgi).*dphi1_g(xgi));
    R_g(2,1) = sum(wi.*dphi1_g(xgi).*dphi0_g(xgi));
    R_g(2,2) = sum(wi.*dphi1_g(xgi).*dphi1_g(xgi));
    
    N = length(xi);   
    M = sparse(N,N);
    R = sparse(N,N);

% ENSAMBLAJE 

    %En los elementos de la diagonal principal se van solapando los valores de las matrices de referencia
    
    for i=1:Ne

     hi=xi(i+1)-xi(i); %Por si hubiera usado nodos no equiespaciados
     M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + hi*0.5*M_g;
     R(i:i+1,i:i+1) = R(i:i+1,i:i+1) + (2/hi)*R_g;
    
    end


