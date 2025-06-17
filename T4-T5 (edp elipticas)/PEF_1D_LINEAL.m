%% PROBLEMA ELÍPTICO 1D ELEMENTOS FINITOS LINEALES. CONDICIONES DIRICHLET HOMOGÉNEAS
% a*u - k*u'' = f
% solo hay que cabiar a,k,f, los nodos (o la h) y el intervalo

%% MODIFICAR

    f = @(x) 0*x+1;
    a=1;
    k=0.1;

% Intervalo
    c=0;
    d=1;

% NODOS:
    % equiespaciados:
    Ne = 25;      %num elementos
    h = (d-c)/Ne; 
    xi = c:h:d;   %nodos
    
    %Para evitar oscilaciones en los extremos se suele utilizar nodos no equiespaciados:
    %xi = [0:1e-3:0.01 0.01:0.05:0.99 0.99:1e-3:1]; %se ponen nodos más juntos en los extremos que por el medio
    %Ne = length(xi)-1; %para el caso lineal
    
%% CREAR MATRICES DE MASA Y RIGIDEZ EN EL ELEMENTO DE REFERENCIA (con cuadratura gauss legendre)

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

%% ENSAMBLAJE 

    %En los elementos de la diagonal principal se van solapando los valores de las matrices de referencia
    
    for i=1:Ne

     hi=xi(i+1)-xi(i); %Por si hubiera usado nodos no equiespaciados
     M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + hi*0.5*M_g;
     R(i:i+1,i:i+1) = R(i:i+1,i:i+1) + (2/hi)*R_g;
    
    end

%% SOLUCIÓN

A = a*M+k*R;
A0 = A; %En vez de redimensionalizar, ponemos cero en la 1 y ult fila y columna (menos en posición 11 y NN)
A0([1 end],:) = 0; %Es lo mismo que poner A0(1,:)=0 y A0(end,:)=0
A0(:,[1 end]) = 0; %Es lo mismo que pober A0(:,1)=0 y A0(:,end)=0
A0(1,1) = 1;
A0(end,end) = 1;

fi=f(xi)';

vect_b=M*fi;
vect_b0=vect_b;
vect_b0([1 end])=0;


uh=A0\vect_b0; % solución del PEF

figure

plot(xi,uh)
