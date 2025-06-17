load("mallado_barra.mat")


xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_d = find( (xi<=0.5 & yi==0) | (xi<=0.5 & yi ==0.2) | (xi>=10.5 & yi==0) | (xi>=10.5 & yi==0.2) | xi==0 | xi==11);

[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R12 M] = assema(p,t,[0 1 0 0]',1,0);
[R21 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);

E = 2.7e8;
nu = 0.2;
rho = 1e3;
g = 9.8;

mu = E*0.5/(1+nu);
lambda = E*nu/((1+nu)*(1-2*mu));

A11 = (lambda+2*mu)*R11 + mu*R22;
A12 = lambda*R21 + mu*R12;
A21 = lambda*R12 + mu*R21;
A22 = (lambda+2*mu)*R22 + mu*R11;

A = [A11 A12; A21 A22];

f1i = 0*xi';
f2i = 0*xi' - rho*g;

vect_b1 = M*f1i;
vect_b2 = M*f2i;

vect_b = [vect_b1; vect_b2];

%condiciones de contorno
fron_d = [fron_d fron_d+length(xi)];
A0 = A;
A0(fron_d,:) = 0;
A0(:,fron_d) = 0;
for i = fron_d
    A0(i,i) = 1;
end
vect_b(fron_d) = 0;

uh = A0\vect_b;
u1h = uh(1:length(xi));
u2h = uh(length(xi)+1:end);

figure(1)
trisurf(elem,xi,yi,u1h)
figure(2)
trisurf(elem,xi,yi,u2h)
figure(3)
trisurf(elem,xi,yi,u2h)
view(2)
shading interp
colorbar
axis equal
figure(4)
trisurf(elem,xi+u1h',yi+u2h',0*xi)%no me sale el plot igual que a el
axis equal