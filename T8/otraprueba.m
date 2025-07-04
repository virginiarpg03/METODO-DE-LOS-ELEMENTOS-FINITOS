% Ecuación del calor
% Trabajo, parte bomba de calor

% Cargar resultados de Stokes desde bomba_calor.m
% bomba_calor

% Parámetros físicos
f = @(x, y, t) 0 * x; % Término fuente nulo
k = 0.025; % Coeficiente de difusión térmica
a = 0; % Coeficiente de reacción, nulo en este caso

% Parámetros temporales
t0 = 0; % Tiempo inicial
tf = 100; % Tiempo final
dt = 0.1; % Paso de tiempo
Nt = (tf - t0) / dt; % Número de pasos de tiempo

% Extraer nodos y elementos
xi = p(1, :); % Coordenadas x de los nodos
yi = p(2, :); % Coordenadas y de los nodos
elem = t(1:3, :)'; % Elementos del mallado

% Identificar fronteras y condiciones de contorno
front_d = unique([e(1, :) e(2, :)]); % Nodos en las fronteras
ind_1 = find(xi == 3.25 & yi >= 0 & yi <= 0.08 | ...
             xi == 4.75 & yi >= 0 & yi <= 0.08 | ...
             yi == 0.08 & xi >= 3.25 & xi <= 4.75); % Radiador
ind_2 = find(yi == 3 & xi >= 0 & xi <= 0.1 | ...
             yi == 4 & xi >= 0 & xi <= 0.1 | ...
             xi == 0.1 & yi >= 3 & yi <= 4); % Ventana
ind_3=find(yi==0 | yi==7 | xi==0 | xi==8);


% Crear matrices de convección
Cx_b1 = sparse(length(xi), length(xi));
Cy_b2 = sparse(length(xi), length(xi));

for i = 1:length(xi)

    % Inicializamos phi_i para el nodo i-ésimo
    phi_i = 0 * xi'; 
    phi_i(i) = 1; % phi_i es 1 en el nodo i y 0 en el resto
    
    % Gradientes de la función phi_i
    [phi_i_x, phi_i_y] = pdecgrad(p, t, 1, phi_i);
    
    % Incorporar b1h (velocidad en x) en la matriz de convección
    b1_local = b1h(i); % Valor de b1h en el nodo i
    [R, M, b] = assema(p, t, 1, 1, phi_i_x * b1_local); % Multiplicamos gradiente por b1h
    Cx_b1(:, i) = b;

    % Incorporar b2h (velocidad en y) en la matriz de convección
    b2_local = b2h(i); % Valor de b2h en el nodo i
    [R, M, b] = assema(p, t, 1, 1, phi_i_y * b2_local); % Multiplicamos gradiente por b2h
    Cy_b2(:, i) = b;

end




% Crear matrices de masa y rigidez
[R, M] = assema(p, t, 1, 1, 0); % Matriz de rigidez y masa

% Ensamblar matriz A para el sistema implícito
A = M + dt * 0.5 * Cx_b1 + dt * 0.5 * Cy_b2 + k * dt * 0.5 * R;

% Ajustar matriz A para condiciones de contorno
A0 = A;
A0(front_d, :) = 0;
A0(:, front_d) = 0;
for i = front_d
    A0(i, i) = 1;
end

% Condiciones iniciales y de contorno
gi = 0 * xi'; % Inicializar con ceros
gi(ind_1) = 5; % Condición de contorno en el radiador
gi(ind_2) = 30; % Condición de contorno en la ventana
gi(ind_3)= 17;  % condiciones paredes

uhn = 0 * xi' +20; % Inicializar temperatura en el dominio

% Visualización inicial
trisurf(elem, xi, yi, uhn);
title('t = 0');
pause;

% Bucle temporal para la evolución de la temperatura
for n = 1:Nt
    % Ensamblar lado derecho del sistema
    vect_b = (M - dt * 0.5 * Cx_b1 - dt * 0.5 * Cy_b2 - k * dt * 0.5 * R) * uhn - A * gi;
    vect_b(front_d) = 0;

    % Resolver sistema lineal
    whn = A0 \ vect_b;
    uhn = whn + gi;

    % Visualizar evolución de la temperatura
    trisurf(elem, xi, yi, uhn);
    title(['t = ', num2str(n * dt)]);
    shading interp
    colorbar
    pause(0.005);
end

% Visualización final
trisurf(elem, xi, yi, uhn);
title('Temperatura final');
