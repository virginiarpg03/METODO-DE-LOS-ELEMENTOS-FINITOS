% interpolar una función mediante elementos triangulares

function [fh,error] = elementos_triangulos(f,V1,V2,V3,X,m)

% f(x,y)= funcion a interpolar
% V1, V2, V3 == vértices del triánglo (en columna)
% X == punto a aproximar (en columna)
% m == grado del polinomio interpolador

X1=V1; X2=V2; X3=V3;
f= f(X(1),X(2)); %valor real

%lineal
    if m==1 

        bi=X1;
        Ai= [X2-X1 , X3-X1];
        Xg= Ai\(X-bi);
        
        xg=Xg(1);
        yg=Xg(2);
        
        % queremos hallar phi(X). Necesitamos la phi gorro ya que ==> phi(X)=phigorro(Xg) 
        
        f1= f(X1(1) , X1(2));
        f2= f(X2(1) , X2(2));
        f3= f(X3(1) , X3(2));
        
        %funciones base
        phi1= 1 -xg -yg;
        phi2= xg;
        phi3= yg;
        
        fh= f1*phi1 + f2*phi2 + f3*phi3; %valor interponado
        error=abs(f-fh);
    end

%cuadrático
    if m==2 
        % también los puntos medios
        X4=0.5*(X1+X2);
        X5=0.5*(X2+X3);
        X6=0.5*(X1+X3);
        
        % para hallar xg, yg es lo mismo que antes
        bi=X1;
        Ai= [X2-X1 , X3-X1];
        Xg= Ai\(X-bi);
        
        xg=Xg(1);
        yg=Xg(2);
        
        % VALORES DE LA FUNCION EN CADA PUNTO
        f1 = f(X1(1) , X1(2));
        f2 = f(X2(1) , X2(2));
        f3 = f(X3(1) , X3(2));
        f4 = f(X4(1) , X4(2));
        f5 = f(X5(1) , X5(2));
        f6 = f(X6(1) , X6(2));
        
        %FUNCIONES BASE
        phi1=+1/4*(xg-1)*(yg-1);      
        phi2=-1/4*(xg+1)*(yg-1); 
        phi3=+1/4*(xg+1)*(yg+1); 
        phi4=-1/4*(xg-1)*(yg+1); 
        
        
        fh=f1*phi1+f2*phi2+f3*phi3+f4*phi4+f5*phi5+f6*phi6;
        error=abs(fh-f);
    end

end
        


