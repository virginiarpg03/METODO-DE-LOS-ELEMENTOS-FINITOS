function [pf] = EF_Recta_Real_Chevy(fun, a0, b0, Ne0, m0)
    f = fun; %Función
    
    a = a0; % Intervalo inferior
    b = b0; % Intervalo superior
    
    Ne = Ne0; % Número de elementos
    
    h = (b-a)/Ne; % Longitud del elemento
     
    m= m0; % Grado del polinomio + 1;
    
    x = a:0.001:b;
    
    xi = zeros(1,15);
    
    for k=1:Ne
        
        for l=1:m
            nodo = cos((2*l+1)/(2*m+2)*pi);
            xi(l) = nodo;
        end
        xi = 5*xi; % Función afín que transforma los nodos
        %xi = linspace(a + h*(k-1), a + h*k, m); % Nodos de cada elemento
        x = a + h*(k-1):0.001:a + h*k;
        n = length(xi);
        pf = 0.*x;
    
        for i = 1:n
        p = ones(1,length(x));
        
        p=p.*f(xi(i));
    
            for j = 1:n
                
                if i ~= j
                    p=p.*(x-xi(j))/(xi(i)-xi(j));
                end
            end
        pf = pf + p;
        end
    
    
    end

end
