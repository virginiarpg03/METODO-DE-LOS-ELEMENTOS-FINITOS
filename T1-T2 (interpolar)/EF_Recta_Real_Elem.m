function [pf] = EF_Recta_Real_Elem(fun, a0, b0, Ne0, m0)
    hold off

    f = fun; %Función
    
    a = a0; % Intervalo inferior
    b = b0; % Intervalo superior
    
    Ne = Ne0; % Número de elementos
    
    h = (b-a)/Ne; % Longitud del elemento
     
    m= m0; % Grado del polinomio + 1;
    
    x = a:0.001:b;
    plot(x, f(x), "k")
   
    hold on
    for k=1:Ne
        
        xi = linspace(a + h*(k-1), a + h*k, m); % Nodos de cada elemento
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
    
        plot(x, pf, "r")
    end

    hold off
end

