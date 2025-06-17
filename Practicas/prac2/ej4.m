%% EJ 4

% hemos hallado uh evaluada en cada nodo (xi,yi), para usar la cuadratura
% de elementos finitos en triángulos los pesos serán w1=w2=w3= 1/6


% Definir puntos de cuadratura en la referencia del triángulo (vértices)
xi=p(1,:);   
yi=p(2,:);    

elem=t(1:3,:)';    %guarda los nodos de cada triángulo

Ne=length(elem);    

Xg=[0 1 0;0 0 1];   
wi=[1/6 1/6 1/6];    

Q=0;

% Bucle sobre cada elemento 
for i = 1:Ne
    % Índices de los nodos del elemento
    nodos = elem(i, :);

    % Valores de u_h en los nodos del elemento
    uh_nodos = uh(nodos);

    X1 = p(:, nodos(1));
    X2 = p(:, nodos(2));
    X3 = p(:, nodos(3));

    Ai = [X2 - X1, X3 - X1];

    %necesitamos multiplicarlo por el área
    detAi = abs(det(Ai));

   % cuadratura
    for j = 1:3
        % Evaluar uh en el punto de cuadratura (coincide con el nodo j)
        f = uh_nodos(j);

        % Acumular la contribución al integral
        Q = Q + wi(j) * detAi * f;
    end
end

Q
