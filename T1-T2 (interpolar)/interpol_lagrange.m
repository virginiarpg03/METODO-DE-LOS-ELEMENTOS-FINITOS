%% función interpolación de lagrange:
%  meto los nodos xi y los valores yi de esos nodos. Si la función es conocida
%  primero hallo los nodos yi de forma:
%  f=@(x)= ... ; yi=f(xi)
%  x: puntos en ls cuales obtener el valor evaluado
%  devuelve el polinomio de interpolacion

function [p] = interpol_lagrange (xi,fi, x)

p = 0*x; 

for i=1:length(xi) 
    li=0*x+1;

    for j=1:length(xi)
        if i ~= j
         li = li.*(x-xi(j))/(xi(i)-xi(j));
        end
    end

    p = p + fi(i)*li;

end