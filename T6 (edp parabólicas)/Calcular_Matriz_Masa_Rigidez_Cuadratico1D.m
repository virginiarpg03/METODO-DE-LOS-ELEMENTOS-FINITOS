%% CREAR MATRICES DE MASA Y RIGIDEZ EN EL ELEMENTO DE REFERENCIA 1D G2 (con cuadratura gauss legendre) Y ENSAMBLARLAS

% Como ahora tenemos polinomios de orden 2, tenemos que tomar las phis de orden 2.
%funciones base y sus derivadas:

    phi0 = @(x) -0.5.*x.*(1-x);
    phi1 = @(x) 1-x.*x;
    phi2 = @(x) 0.5.*x.*(1+x);
    
    dphi0 = @(x) -0.5+x;
    dphi1 = @(x) -2*x;
    dphi2 = @(x) 0.5+x;


% cuadratura de Gauss Legendre de 3 puntos:
xgi = [-sqrt(3/5) 0 sqrt(3/5)];
wi = [5/9 8/9 5/9];

%definimos las matrices de referencia, que ahora son de 3x3
    Mg = zeros(3,3);
    
    phi0i = phi0(xgi);
    phi1i = phi1(xgi);
    phi2i = phi2(xgi);
    
    Mg(1,1) = sum(wi.*phi0i.*phi0i);
    Mg(2,2) = sum(wi.*phi1i.*phi1i);
    Mg(3,3) = sum(wi.*phi2i.*phi2i);
    Mg(1,2) = sum(wi.*phi1i.*phi0i);
    Mg(2,1) = Mg(1,2);                  %M es simétrica
    Mg(1,3) = sum(wi.*phi0i.*phi2i);
    Mg(3,1) = Mg(1,3);
    
    Mg(2,3) = sum(wi.*phi1i.*phi2i);
    Mg(3,2) = Mg(2,3);


%matriz Rg
    Rg=zeros(3,3);
    
    dphi0i = dphi0(xgi);
    dphi1i = dphi1(xgi);
    dphi2i = dphi2(xgi);
    
    Rg(1,1) = sum(wi.*dphi0i.*dphi0i);
    Rg(2,2) = sum(wi.*dphi1i.*dphi1i);
    Rg(3,3) = sum(wi.*dphi2i.*dphi2i);
    Rg(1,2) = sum(wi.*dphi1i.*dphi0i);
    Rg(2,1) = Rg(1,2);                  %simétrica
    Rg(1,3) = sum(wi.*dphi0i.*dphi2i);
    Rg(3,1) = Rg(1,3);
    
    Rg(2,3) = sum(wi.*dphi1i.*dphi2i);
    Rg(3,2) = Rg(2,3);

%calculamos las matrices M y R
    M = sparse(2*Ne+1, 2*Ne+1);
    R = sparse(2*Ne+1, 2*Ne+1);
    
    for i=1:Ne
        hi = xi(i+1) - xi(i); %En este caso vale h
        j = 2*(i-1)+1;
        M(j:j+2, j:j+2) = M(j:j+2, j:j+2) + hi*0.5*Mg;
        R(j:j+2, j:j+2) = R(j:j+2, j:j+2) + 2*Rg/hi;
    end
