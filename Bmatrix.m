function B=Bmatrix(n,j,l,r,k)

ept=0.0001;
B1=@(x)4*k*1i*((((x.^2-1/2).^2)./(x.*(1-x.^2).^(1/2))).*shankel(n,j,k*x)).*sbesselj(n,l,k*x);
B2=@(x)1/(2*n+4*j+3)*(delta(j,l))./((x.^2).*(1-x.^2).^(1/2));
B3=@(x)4*k*1i*((x.*(r^2-x.^2).^(1/2)).*shankel(n,j,k*x)).*sbesselj(n,l,k*x);

B=integral(@(x)B1(x),ept,1)+integral(@(x)B2(x),ept,1)+integral(@(x)B3(x),ept,r);