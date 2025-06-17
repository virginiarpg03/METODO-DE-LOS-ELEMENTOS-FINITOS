
f=@(x) (3+x).*cos(pi*x*0.25).^2 ;

xi= [0 1 3];
fi=f(xi);

x=[2 2.4 3.5 4];

p= interpol_lagrange(xi ,fi,x)

