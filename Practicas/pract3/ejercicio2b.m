%% EJERCICIO 2 b

% la cuadratura de simpsom : Qf= h/3*(f(a) + 4*f(a+h) + f(b)),

% sacar los nodos de las fronteras 3 y 4 y valores de uh en esos nodos
%save("mallado_stokes","t","e","p")

load("")
xi=p(1,:);   
yi=p(2,:);  

% frontera 3 y 4
front3= find(yi==5 | xi>=8.5);
front4= find(xi==15 | yi>=-1.5);

% Extraer los valores de u en las fronteras
u3 = uh(front3); 
u4 = uh(front4); 

% voy cogiendo de 2 en 2 nodos
Q=0
for i=1: (length(front3)-1)
    a = front3(i);
    b = front3(i+1);
    h = (b-a)*0.5;
    ua= u3(i);
    ub= u3(i+1);
    um= (ub-ua)*0.5;

    % cuadratura:
    Q= Q * h/3*(ua + 4*um + ub);

end

Q3= Q



Q=0
for i=1: (length(front3)-1)
    a = front4(i);
    b = front4(i+1);
    h = (b-a)*0.5;
    ua= u4(i);
    ub= u4(i+1);
    um= (ub-ua)*0.5;

    % cuadratura:
    Q= Q * h/3*(ua + 4*um + ub);

end

Q4= Q









