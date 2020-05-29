function A=Amatrix(m,j,l,r,k)
 
 ept=0.0001;
 A1=@(x)4*k*1i*((((x.^2-1/2).^2)./(x.*(r^2-x.^2).^(1/2))).*shankel(m,j,k*x)).*sbesselj(m,l,k*x);
 A2=@(x)1/(2*m+4*j+3)*(delta(j,l))./((x.^2).*(r^2-x.^2).^(1/2));
 A3=@(x)4*k*1i*((x.*(1-x.^2).^(1/2)).*shankel(m,j,k*x)).*sbesselj(m,l,k*x);
        
 A=integral(@(x)A1(x),ept,r)+integral(@(x)A2(x),ept,r)+integral(@(x)A3(x),ept,1);