function [] = Acota_Errores(a0, b0, n0)

    n = n0; % Número de nodos/Grado de la inter-1
    
    a = a0; %Extremo inicial
    
    b = b0; %Extremo final
    
    f = @(x) log(x); % Función a interpolar
    
    df0 = @(x) 2*x.^(-3);
    dfneg0 = @(x) -2*x.^(-3); %Derivadas de grado n
    
    w0 = @(x) (x-9).*(x-9.5).*(x-10); % Multiplicación de Nodos
    w0neg = @(x) -(x-9).*(x-9.5).*(x-10);
    
    [r, minimow0] = fminbnd(w0, a, b);
    [r, minimow0neg] = fminbnd(w0neg, a, b);
    
    
    w = max(abs(minimow0), abs(minimow0neg)); % Máximo de la funcion w
    
    [r, minimodf] = fminbnd(df0, a, b);
    [r, minimodfneg] = fminbnd(dfneg0, a, b);
    
    d3f = max(abs(minimodf), abs(minimodfneg)); %Máximo de la derivada n-esima
    
    
    E = d3f*w/factorial(n); % Error máximo
    
    fprintf("El error máximo es: %6f", E)

end