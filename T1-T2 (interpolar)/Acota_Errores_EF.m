function [] = Acota_Errores_EF(n0, h)

    n = n0; % Número de nodos/Grado de la inter+1
    
    a = -1; %Extremo inicial
    
    b = 1; %Extremo final
    
    f = @(x) log(x); % Función a interpolar
    
    df0 = @(x) 1;
    dfneg0 = @(x) -1; %Derivadas de grado n
    
    w0 = @(x) (x+1).*(x).*(x-1); % Multiplicación de Nodos
    w0neg = @(x) -(x+1).*(x).*(x-1);
    
    [r, minimow0] = fminbnd(w0, a, b);
    [r, minimow0neg] = fminbnd(w0neg, a, b);
    
    
    w = max(abs(minimow0), abs(minimow0neg)); % Máximo de la funcion w
    
    [r, minimodf] = fminbnd(df0, a, b);
    [r, minimodfneg] = fminbnd(dfneg0, a, b);
    
    df = max(abs(minimodf), abs(minimodfneg)); %Máximo de la derivada n-esima
    
    
    E = h^n*df*w/factorial(n)/2^n; % Error máximo
    
    fprintf("El error máximo es: %6f", E)

end