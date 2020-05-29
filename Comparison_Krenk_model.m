clear;
close;
clc;

Km=26;%bulk modulus of the background medium,GPa
Gm=31;%shear modulus of the background medium
por=0.1;%porosity of the background medium
Ks=37;%bulk modulus of the solid 
Kf=0.00001;%bulk modulus of the fluid
dens=2.65;%density of the solid
denf=0.001;%density of the fluid
vis=0.001*10^(-9);%viscosity of fluid, Gpa*s
perm=10^(-16);%permeability of background, m^2
d=1;%diameter of the fracture,m
r=d/2;%radius of the fracture,m
cden=0.05;%fracture density

ep=0.0001;
h=0.001;
Ang=[ep pi/4 pi/2-ep];
Da=length(Ang);
f1=2:0.2:3;
f2=3.02:0.02:4;
f=[f1 f2];
f=10.^f;
Df=length(f);
v=zeros(Df,Da);
inQ=zeros(Df,Da);

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;
den=dens*(1-por)+denf*por;

for AN=1:Da
    
    ang=Ang(AN);

for FN=1:Df

FN

ft=f(FN);
w=2*pi*ft;
 
abk=(2*pi)/d*5;

if(ft<=1000)
DM=10000;
else
DM=5000;
end

k=abk/DM:abk/DM:abk;

wb=por*vis/(denf*perm*10^(-6));
permb=perm*(sqrt(1-1i*w/(2*wb))-1i*w/wb)^(-1);
denb=(1i*vis)/(w*permb)*10^6;
b=(den*M+denb*H-2*denf*C)/(2*(M*H-C^2));
%s1=sqrt(b-sqrt(b^2-(den*denb-denf^2)/(M*H-C^2)));
s1=sqrt((den*denb-denf^2)/(M*H-C^2)/(b+sqrt(b^2-(den*denb-denf^2)/(M*H-C^2))));
s2=sqrt(b+sqrt(b^2-(den*denb-denf^2)/(M*H-C^2)));
s3=sqrt((den*denb-denf^2)/(Gm*denb));

k1=w*s1*10^(-3);
k2=w*s2*10^(-3);
k3=w*s3*10^(-3);

x1=-(H*s1^2-den)/(C*s1^2-denf);
x2=-(H*s2^2-den)/(C*s2^2-denf);
x3=-denf/denb;

kN=length(k);

et1=zeros(1,kN);
et2=zeros(1,kN);
et3=zeros(1,kN);

for lk=1:kN
    
kk=k(lk);
if kk<=k1
    et1(lk)=-1i*sqrt(k1^2-kk^2);
else
    et1(lk)=sqrt(kk^2-k1^2);
end
    et2(lk)=-1i*sqrt(k2^2-kk^2);
if kk<=k3
    et3(lk)=-1i*sqrt(k3^2-kk^2);
else
    et3(lk)=sqrt(kk^2-k3^2);
end

end

% normal fracture discontinuity

E1=-1*(x3-x2)*(H-Gm-C+C*x1-M*x1)*k1^2;
E2=(x3-x1)*(H-Gm-C+C*x2-M*x2)*k2^2;
E=2*(E1+E2);

F1=2*Gm*k.^2-(H-C-M*x1+C*x1)*k1^2-(4*Gm*(et1.*et3).*k.^2)./(2*k.^2-k3^2);
F2=(2*Gm*k.^2-(H-C-M*x2+C*x2)*k2^2-(4*Gm*(et2.*et3).*k.^2)./(2*k.^2-k3^2)).*...
    (((2*(x3-x1)*k.^2+(1+x1)*k3^2).*et1+h/(2*Kf)*(C+M*x1)*k1^2*(2*k.^2-k3^2))./((2*(x3-x2)*k.^2+(1+x2)*k3^2).*et2+h/(2*Kf)*(C+M*x2)*k2^2*(2*k.^2-k3^2)));
Hf=E^(-1)*((((F1-F2).*(2*(x3-x2)*k.^2+(1+x2)*k3^2)).*et1.^(-1)).*k.^(-1))-1;

N=50;
ds=r/N;
Ma=zeros(N,N);
V=zeros(N,1);

Mx=0;
sumM=0;

for mf=0:Mx

if(mf==0)
    Em=1;
else
    Em=2;
end

for I=1:N
    for J=1:N
       si=(I-1)*ds;
       sj=(J-1)*ds; 
       
       HK=((k.*Hf).*besselj(mf+0.5,si*k)).*besselj(mf+0.5,sj*k);
       My=trapz(k,HK)*(si*sj)^0.5;
       if(I==J)
           Ma(I,J)=My*ds+1;
       else
           Ma(I,J)=My*ds;
       end
       
    end
    
    D1=(H-C+C*x1-M*x1)/Gm;
    
    if(si==0)
      if(mf>0)
          Tmv=0;
      else
          Tmv=sqrt(2/pi);
      end
    else
       Tmv=besselj(mf+0.5,k1*si*sin(ang))/sqrt(k1*si*sin(ang)); 
    end
    
    p0=Gm*k1*(D1-2*sin(ang)^2)*Em*(1i)^(mf+1)*sqrt(pi/2)*Tmv;
    V(I)=-p0*si;
end

P=Ma\V;
sum=0;
for I=1:N
    si=(I-1)*ds;
    sum=sum+si^0.5*P(I)*besselj(mf+0.5,si*k1*sin(ang))*ds;
end

sum=sum*sqrt(2*k1*sin(ang)/pi);
sumM=sumM+sum*(-1i)^(mf+1);

end

f0=-(2*k1^2*sin(ang)^2*(x3-x2)+k3^2*(1+x2))/(sin(ang)*E)*sumM;

% shear fracture discontinuity

E=(C^2-(H-Gm)*M)*(x1-x2)*k1^2/(Gm*(C+M*x2));

F1=2*et1-(2*k.^2-k3^2)./et3+(k1^2*(2*k.^2-k3^2)*(H+C*x1))./((2*Gm*et3).*(k.^2));
F2=(k1^2*(C+M*x1))/(k2^2*(C+M*x2))*(2*et2-(2*k.^2-k3^2)./et3+(k2^2*(2*k.^2-k3^2)*(H+C*x2))./((2*Gm*et3).*(k.^2)));
Hf=(F1-F2).*k/E-1;

N=50;
ds=r/N;
Ma=zeros(N,N);
V=zeros(N,1);

Mx=0;
sumM=0;

for mf=0:Mx

if(mf==0)
    Em=1;
else
    Em=2;
end

for I=1:N
    for J=1:N
       si=(I-1)*ds;
       sj=(J-1)*ds; 
       
       HK=((k.*Hf).*besselj(mf+0.5,si*k)).*besselj(mf+0.5,sj*k);
       My=trapz(k,HK)*(si*sj)^0.5;
       if(I==J)
           Ma(I,J)=My*ds+1;
       else
           Ma(I,J)=My*ds;
       end
       
    end
    
    if(si==0)
      if(mf>0)
          Tmv=0;
      else
          Tmv=sqrt(2/pi);
      end
    else
       Tmv=besselj(mf+0.5,k1*si*sin(ang))/sqrt(k1*si*sin(ang)); 
    end
    
    p0=2*cos(ang)*Em*(1i)^mf*sqrt(pi/2)*Tmv;
    V(I)=p0*si;
end

P=Ma\V;
sum=0;
for I=1:N
    si=(I-1)*ds;
    sum=sum+si^0.5*P(I)*besselj(mf+0.5,si*k1*sin(ang))*ds;
end

sum=sum*sqrt(2*k1*sin(ang)/pi);
sumM=sumM+sum*(-1i)^mf;

end

f1=(k1^3*sin(2*ang))/(2*E)*sumM;

fT=f0+f1;

n0=cden/r^3;
keff=k1*(1+(2*pi*n0)/k1^2*fT);
v(FN,AN)=(2*pi*ft)/real(keff);
inQ(FN,AN)=2*imag(keff)/real(keff);

end

end

%Krenk-Schmidt method

vr=(3*Km-2*Gm)/(2*(3*Km+Gm));%possion ratio of the background medium

epa=0.0001;
Ang=[0+epa,pi/4,pi/2-epa];
An=length(Ang);

N=15;
ad=0.95;
suma=(1-ad^N)/(1-ad);
ft=zeros(1,N+1);
ft0=2;
ft1=4;
ftd=ft1-ft0;
t=ftd/suma;

tft=0;
ft(1)=ft0;

for tf=1:N
    tft=tft+ad^(tf-1)*t;
    ft(tf+1)=ft0+tft;
end

ft=10.^ft;
fn=length(ft);
vp=zeros(fn,An);
inQ1=zeros(fn,An);

for AN=1:An
    
  ang=Ang(AN);

for FN=1:fn

FN
f=ft(FN);
w=2*pi*f;
h=sqrt((1-2*vr)*w^2*r^2*den/(2*(1-vr)*Gm))*10^(-3);
k=sqrt(w^2*r^2*den/Gm)*10^(-3);
rb=sqrt(1/2*(1-2*vr)/(1-vr));
LN=5;
JN=5;

mx1=0;
sum0=0;

for m=0:1:mx1

if(m==0)
    em=1;
else
    em=2;
end

A=zeros(JN,LN);
S=zeros(1,LN);

for J=1:JN
    for L=J:LN
        j=J-1;
        l=L-1;
        
        A(J,L)=Amatrix(m,j,l,rb,k);
        A(J,L)=(-1)^(j+l)*A(J,L);
        A(L,J)=A(J,L);
        
    end
end

for L=1:LN
    
    l=L-1;
    S(L)=1/sin(ang)*(rb^(-2)-2*sin(ang)^2)*em*1i^(m+1)*(-1)^l*sbesselj(m,l,h*sin(ang));
    
end

W=S/A;

sum1=0;
s=h*sin(ang);
for J=1:JN
    j=J-1;
    sum1=sum1+(-1)^j*W(J)*sbesselj(m,j,s);
end

sum1=2/k^2*sum1;

if(s>h)
    as=(s^2-h^2)^(1/2);
else
    as=1i*(h^2-s^2)^(1/2);
end

Am1=(s^2-1/2*k^2)*sum1/(s*as);
sum0=sum0+(-1i)^m*cos(ang)*Am1;
end

f0=h^2*r*sum0;

mx2=1;
sum2=0;
for m=0:1:mx2

if(m==0)
    em=1;
else
    em=2;
end

if(m==0)
B=zeros(JN,LN);
T=zeros(1,LN);
n=1;

for J=1:JN
    for L=J:LN
        j=J-1;
        l=L-1;
        
        B(J,L)=Bmatrix(n,j,l,rb,k);
        B(J,L)=(-1)^(j+l)*B(J,L);
        B(L,J)=B(J,L);
        
    end
end

for L=1:LN
    
    l=L-1;
    T(L)=2*(-1)^(l+1)*cos(ang)*sbesselj(1,l,h*sin(ang));
    
end

U=T/B;

sum3=0;
s=h*sin(ang);
for J=1:JN
    j=J-1;
    sum3=sum3+(-1)^j*2*U(J)*sbesselj(1,j,s);
end

Am2=sum3/k^2;

else
    
n1=m+1;
n2=m-1;
B1=zeros(JN,LN);
C1=zeros(JN,LN);
B2=zeros(JN,LN);
C2=zeros(JN,LN);
B3=zeros(JN,LN);
C3=zeros(JN,LN);
B4=zeros(JN,LN);
C4=zeros(JN,LN);

for J=1:JN
    for L=J:LN
        j=J-1;
        l=L-1;
        tp=(-1)^(j+l);
        
        B1(J,L)=Bmatrix(n1,j,l,rb,k);
        B1(J,L)=tp*B1(J,L);
        B1(L,J)=B1(J,L);
        
        C1(J,L)=Cmatrix(n1,j,l,k);
        C1(J,L)=tp*C1(J,L);
        C1(L,J)=C1(J,L);
        
        B3(J,L)=Bmatrix(n2,j,l,rb,k);
        B3(J,L)=tp*B3(J,L);
        B3(L,J)=B3(J,L);
        
        C3(J,L)=Cmatrix(n2,j,l,k);
        C3(J,L)=tp*C3(J,L);
        C3(L,J)=C3(J,L);
        
     end
end

for J=1:JN
    for L=1:LN
        j=J-1;
        l=L-1;
        tp=(-1)^(j+l);
        
        if((j-1)<=l)
            B2(J,L)=tp*Bmatrix(n1,j-1,l,rb,k);
            C2(J,L)=tp*Cmatrix(n1,j-1,l,k);
        else
            B2(J,L)=tp*Bmatrix(n1,l,j-1,rb,k);
            C2(J,L)=tp*Cmatrix(n1,l,j-1,k);
        end
        
        if((j+1)<=l)
            B4(J,L)=tp*Bmatrix(n2,j+1,l,rb,k);
            C4(J,L)=tp*Cmatrix(n2,j+1,l,k);
        else
            B4(J,L)=tp*Bmatrix(n2,l,j+1,rb,k);
            C4(J,L)=tp*Cmatrix(n2,l,j+1,k);
        end
        
    end

end

BM=[B1+C1,-(B4-C4);-(B2-C2),B3+C3];
BM=1/2*BM;

T1=zeros(1,LN);
T2=zeros(1,LN);

for L=1:LN
    l=L-1;
    
    T1(L)=2*(-1)^(l+1)*cos(ang)*em*(1i)^m*sbesselj(n1,l,h*sin(ang));
    T2(L)=-2*(-1)^(l+1)*cos(ang)*em*(1i)^m*sbesselj(n2,l,h*sin(ang));
end

TM=[T1,T2];

UM=TM/BM;

U1=UM(1:JN);
U2=UM(JN+1:2*JN);

sum4=0;
s=h*sin(ang);
for J=1:JN
    j=J-1;
    
    sum4=sum4+(-1)^j*(U1(J)*sbesselj(n1,j,s)-U2(J)*sbesselj(n2,j,s));
end

Am2=sum4/k^2;

end

sum2=sum2+(-1i)^m*cos(ang)*Am2;

end

f1=h^2*r*sum2;

fT=-(f0+f1);

n0=cden/r^3;
k1=h/r;
keff=k1*(1+(2*pi*n0)/k1^2*fT);
vp(FN,AN)=(2*pi*f)/real(keff);
inQ1(FN,AN)=-2*imag(keff)/real(keff);

end

end

f1=2:0.2:3;
f2=3.02:0.02:4;
f=[f1 f2];
f=10.^f;

figure(1)
semilogx(f,v(:,1),'m','linewidth',2);
hold on;
semilogx(ft,vp(:,1),'mo','MarkerSize',8);
hold on;
semilogx(f,v(:,2),'r','linewidth',2);
hold on;
semilogx(ft,vp(:,2),'ro','MarkerSize',8);
hold on;
semilogx(f,v(:,3),'g','linewidth',2);
hold on;
semilogx(ft,vp(:,3),'go','MarkerSize',8);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('Vp(0^{\circ})','Vp(0^{\circ})(K-S)','Vp(45^{\circ})','Vp(45^{\circ})(K-S)','Vp(90^{\circ})','Vp(90^{\circ})(K-S)');
set(gca, 'FontSize', 12);
xlim([100 10000]);
set(gca,'XTick',[100 1000 10000]);

figure(2)
loglog(f,inQ(:,1),'m','linewidth',2);
hold on;
loglog(ft,inQ1(:,1),'mo','MarkerSize',8);
hold on;
loglog(f,inQ(:,2),'r','linewidth',2);
hold on;
loglog(ft,inQ1(:,2),'ro','MarkerSize',8);
hold on;
loglog(f,inQ(:,3),'g','linewidth',2);
hold on;
loglog(ft,inQ1(:,3),'go','MarkerSize',8);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Q_{p}^{-1}','FontSize',12);
legend('Vp(0^{\circ})','Vp(0^{\circ})(K-S)','Vp(45^{\circ})','Vp(45^{\circ})(K-S)','Vp(90^{\circ})','Vp(90^{\circ})(K-S)');
set(gca, 'FontSize', 12);
xlim([100 10000]);
set(gca,'XTick',[100 1000 10000]);








