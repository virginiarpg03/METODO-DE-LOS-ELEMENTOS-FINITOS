function [p] = elementos_finitos_1D (f,a,b,Ne,m)
% f== funcion a aproximar para obtener los fi; si me dan fi pues improviso
% (a,b) == intervalo que voy a interpolar
% Ne == número elementos finitos                Ne=(b-a)/h
% m == grado del polinomio interpolador

h = (b-a)/Ne; %distancia entre elementos

figure ('Name', 'aproximación elementos finitos')

for i=1:Ne
    c=a+h;
    xi=linspace(a,c,m+1);
    fi=f(xi);
    x=a:0.001:c;
    p=interpol_lagrange(xi,fi,x);

    hold on
    plot (x,f(x),'k') %real
    plot (x,p)
    hold off
    a=c; %avanzar
end













