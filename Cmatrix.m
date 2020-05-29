function C=Cmatrix(n,j,l,k)

ept=0.0001;
C1=@(x)k*1i*(((x.^(-1)).*(1-x.^2).^(1/2)).*shankel(n,j,k*x)).*sbesselj(n,l,k*x);
C2=@(x)(1/(2*n+4*j+3)*delta(j,l))./((x.^2).*(1-x.^2).^(1/2));

C=integral(@(x)C1(x),ept,1)+integral(@(x)C2(x),ept,1);

