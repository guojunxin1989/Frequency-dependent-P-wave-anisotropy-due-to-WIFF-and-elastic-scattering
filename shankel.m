function h2=shankel(m,j,P)
N=m+2*j+1;
j=(sqrt(pi./(2*P))).*besselj(N+1/2,P);
y=(-1)^(N+1)*(sqrt(pi./(2*P))).*besselj(-N-1/2,P);
h2=j-1i*y;