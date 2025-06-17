clear all
close all
clc 

f = @(x) x.^3;

a=0;
b=10;

%% TRAPECIOS (PRECISIÓN 1)
h= b-a;
Qf= h*0.5*(f(b)+f(a))

%% TRAPECIOS COMPUESTA (par o impar)
clear n h Qf
n=10; %número de particiones
h= (b-a)/n ;
sumf= 0;

for i= 1:n-1
    x(i)= a+h*i;
    sumf= sumf + f(x(i));
end

Qf= h*0.5*(f(a) + 2*sumf + f(b))

%% REGRA DE SIMPSON (PRECISIÓN 3)
clear n h Qf
h = (b-a)*0.5;
Qf= h/3*(f(a) + 4*f(a+h) + f(b))

% %% REGLA 3/8 DE SIMPSON
% clear n h Qf
% h = (b-a)*0.5;
% Qf= h*3/8*(f(a)+3*f(a+h)+3*f(a+2*h)+f(b))

%% REGLA SIMPSON COMPUESTA                             ALGO MAL?!!!!!
clear n h Qf
n=2; %número de particiones
h = (b-a)/n;
s1=0;
s2=0;
for i=1:n-1
    x=a +i*h;

    if rem(i,2)==0
       s2 = s2 + f(x);
    else
       s1 = s1 + f(x);
    end
end

Qf = h/3 * (f(a) + 4*s1 + 2*s2 + f(b))


