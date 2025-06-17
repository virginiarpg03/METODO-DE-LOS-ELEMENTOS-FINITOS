function [val_interp,val_real, error]=elementos_cuadrados(f,V1,V2,V3,V4,Xint,m)     

%f es la funcion a interpolar 
%(V1,V2,V3,V4) son los vertices del rectangulo 
%Xint es el punto donde vamos a interpolar 
%m es el grado de los polinomios de interpolacion 

val_real=f(Xint(1),Xint(2)); 
    if m==1        
        Ai=[V2(1)-V1(1) 0; 0 V3(2)-V2(2)];  %Matriz de la transformacion 
        bi=[V2(1)+V1(1); V3(2)+V2(2)]; 
        Xgorro=Ai\(2*Xint-bi); 
        xg=Xgorro(1); yg=Xgorro(2); 
         
        phi1=+1/4*(xg-1)*(yg-1);      
        phi2=-1/4*(xg+1)*(yg-1); 
        phi3=+1/4*(xg+1)*(yg+1); 
        phi4=-1/4*(xg-1)*(yg+1); 
        
        f1=f(V1(1),V1(2)); 
        f2=f(V2(1),V2(2)); 
        f3=f(V3(1),V3(2)); 
        f4=f(V4(1),V4(2));  
        val_interp=phi1*f1+phi2*f2+phi3*f3+phi4*f4; 
        error=abs(val_interp-val_real);

    end 
    if m==2        
        Ai=[V2(1)-V1(1) 0; 0 V3(2)-V2(2)];  
        bi=[V2(1)+V1(1); V3(2)+V2(2)]; 
        Xgorro=Ai\(2*Xint-bi); 
        xg=Xgorro(1); yg=Xgorro(2); 
         
        V5=(V1+V2)/2; V6=(V2+V3)/2; 
        V7=(V3+V4)/2; V8=(V4+V1)/2; 
        V9=(V7+V5)/2;

        phi1=0.25*(xg^2-xg)*(yg^2-yg);      
        phi2=0.25*(xg^2+xg)*(yg^2-yg); 
        phi3=0.25*(xg^2+xg)*(yg^2+yg); 
        phi4=0.25*(xg^2-xg)*(yg^2+yg); 
        phi5=-0.5*(xg^2 -1)*(yg^2-yg); 
        phi6=-0.5*(xg^2+xg)*(yg^2 -1); 
        phi7=-0.5*(xg^2 -1)*(yg^2+yg); 
        phi8=-0.5*(xg^2-xg)*(yg^2 -1); 
        phi9=(xg^2-1)*(yg^2-1); 
        
        f1=f(V1(1),V1(2)); f2=f(V2(1),V2(2)); f3=f(V3(1),V3(2)); %Valor en los vertices 
        f4=f(V4(1),V4(2)); f5=f(V5(1),V5(2)); f6=f(V6(1),V6(2));  
        f7=f(V7(1),V7(2)); f8=f(V8(1),V8(2)); f9=f(V9(1),V9(2)); 
        
        val_interp=phi1*f1+phi2*f2+phi3*f3+phi4*f4+phi5*f5+phi6*f6+phi7*f7+phi8*f8+phi9*f9; 
        error=abs(val_interp-val_real);
    end 
end

