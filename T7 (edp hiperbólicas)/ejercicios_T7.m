%% EJERCICIOS TEMA 7: EDPS HIPERBÓLICAS


%---------- PROBLEMA HIPERBÓLICO -----------%
%
% utt - 2*k*ut - c2*lapl(u) = f
%
% u(gamma_d) = 0
% u(x,0) = u0
% ut(x,0) = v0
%
%------- discretización: NEWMARK --------%
%
% suponemos k=0, beta!=0
%
% a_0 = c2*lapl(u_0) + f_0
% u_n+1 = u_n + dt*v_n + 0.5*dt^2*(2*beta*a_n+1 + (1-2beta)*a_n)
% a_n+1 se despeja de lo anterior
% v_n+1 = v_n + dt*(gamm*a_n+1 + (1-gamma)*a_n)
%
%------- forma matricial --------%
%
% A = M  + (beta*c2*dt^2)R
%
% M*a_0 = -R*u_0 + M*f_0
% A0*u_n+1 = M(u_n + dt*v_n + 0.5*dt^2*(1-2*beta)*a_n + beta*dt^2*f_n+1)
% beta*dt^2 * a_n+1 = u_n+1 - u_n - dt*v_n - 0.5*dt^2*(1-2*beta)*a_n)
% v_n+1 = v_n + dt*(gamma*a_n+1 + (1-gamma)*a_n)


%% EJERCICIO 7

% utt + 2kut -c^2*uxx = f
% u(x,0) = h0(x)
% ut(x,0) = j0(x)

load("malla_ej6.mat")

% Parámetros problema
c2 = 1;  
% f = @(x,y,t) 0*x;
k = 0.005;
u0 = @(x,y) 0*x;
v0 = @(x,y) x*0;

% Parámetros cuadraturas

             % trapecios             (implícito) beta = 0.25 , gamma= 0.5)
             % diferencias centradas (explícito) beta = 0    , gamma= 0.5)

gamma = 0.5;  
beta = 0.25; 

% mallado 
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_tot = unique([e(1,:) e(2,:)]);

% Partición temporal 
t0=0;
tf=2;

Nt = 1000;
dt = (tf-t0)/Nt;

%  matrices e inicializo mis condiciones iniciales
[R M]=assema(p,t,1,1,0);

 %fi = f(xi,yi,0)';
 ui = u0(xi,yi)';
 vi = v0(xi,yi)';  
 
% en este ejercicio la funcion está definida a trozos
 f0=zeros(1, length(xi)); 
for i=1:length(xi)
    f0(i)=fun_trozos(xi(i),yi(i),0);
end 
f0=f0';

%PINTO LA U INCIAL
trisurf(t(1:3,:)',xi,yi,ui)
xlabel('x'); 
ylabel('y'); 
zlabel('u(x,y,t)');
title('posición incial');
pause

% PINTO V INICIAL
trisurf(elem,xi,yi,vi)
xlabel('x'); 
ylabel('y'); 
zlabel('v(x,y,t)');
title('velocidad incial');
pause

    % PASO 1: OBTENER a0
    A = M;
    vect_b = M*f0-c2*R*ui;

    ai = A\vect_b;


    % RESOLVER PROBLEMAS ELÍPTICOS PARA CADA t 
    
    A = M + c2*beta*dt^2*R; 
    
    A0 = A;
    A0(fron_tot,:)=0;
    A0(:,fron_tot)=0;
    for i=fron_tot
        A0(i,i)=1;
    end    


    for i=1:Nt

        % Elíptico principal

        ui_ant = ui;
        vi_ant = vi;
        ai_ant = ai;
        
        fi=zeros(1, length(xi)); 

        for i=1:length(xi)
            fi(i)=fun_trozos(xi(i),yi(i),(n-1)*dt);
        end 
        fi=fi';
        
        vect_b = M*ui_ant + dt*M*vi_ant +  dt^2*0.5*(1-2*beta)*M*ai_ant  + beta*dt^2*M*fi;
        vect_b0=vect_b;
        vect_b0(fron_tot)=0;
    
        ui = A0\vect_b0; 
        ai = (1/(beta*dt^2)) * (ui - ui_ant - dt*vi_ant - dt^2*0.5*(1-2*beta)*ai_ant); 
        vi = vi_ant + dt*(gamma*ai + (1-gamma)*ai_ant); 

        
        % Representacion gráfica

        trisurf(t(1:3,:)',xi,yi,ui )
        title(n*dt);
        xlabel('x'); 
        ylabel('y'); 
        zlabel('u(x,y,t)');
        pause
  end

%%
function[f]=fun_trozos(x,y,t)
    if (x.^2 + y.^2)< 0.01
            f=@(x,y,t) sin(2*pi*t); 
            f=f(x,y,t);
    else 
            f=0; 
    end 
end

%% EJERCICIO 4