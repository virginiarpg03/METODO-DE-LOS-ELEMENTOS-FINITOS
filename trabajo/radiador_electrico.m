 %Trabajo, parte radiador eléctrico

f=@(x,y,t) 0*x;

a=0;   
k=0.025;   %coef de difusión térmica del aire para una temperatura ambiente de 20 grados

t0=0;
tf=100;
dt=0.1;
Nt=(tf-t0)/dt;

xi=p(1,:);
yi=p(2,:);
elem=t(1:3,:)';

front_d=unique([e(1,:) e(2,:)]);
ind_1=find(xi==3.25 & yi>=0 & yi<=0.08 | xi==4.75 & yi>=0 & yi<=0.08 | yi==0.08 & xi>=3.25 & xi<=4.75);
ind_2=find(yi==3 & xi>=0 & xi<=0.1 | yi==4 & xi>=0 & xi<=0.1 | xi==0.1 & yi>=3 & yi<=4);
ind_3=find(yi==0 | yi==7 | xi==0 | xi==8);


[R M]=assema(p,t,1,1,0);  
A= M + k*dt*0.5*R;
A0=A;
A0(front_d,:)=0;
A0(:,front_d)=0;

for i=front_d
    A0(i,i)=1;
end


gi=0*xi';
gi(ind_1)=5;
gi(ind_2)=30; 
gi(ind_3)=17;

uhn=0*xi'+20;
trisurf(elem,xi,yi,uhn)
title(0)
pause
for n=1:Nt
    vect_b= M*uhn - k*dt*0.5*R*uhn - A*gi;
    vect_b(front_d)=0;

    whn=A0\vect_b;
    uhn=whn+gi;
    trisurf(elem,xi,yi,uhn)
    title(n*dt)
    view(2)
    shading interp
    colorbar
    pause (0.005)
end