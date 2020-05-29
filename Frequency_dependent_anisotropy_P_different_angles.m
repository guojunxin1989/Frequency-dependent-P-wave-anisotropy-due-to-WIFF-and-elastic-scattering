clear;
close;
clc;

Km=26;%bulk modulus of the background medium,GPa
Gm=31;%shear modulus of the background medium
por=0.1;%porosity of the background medium
Ks=37;%bulk modulus of the solid 
Kf=2.25;%bulk modulus of the fluid
dens=2.65;%density of the solid
denf=1;%density of the fluid
vis=0.001*10^(-9);%viscosity of fluid, Gpa*s
perm=10^(-15);%permeability of background, m^2
d=1;%diameter of the fracture,m
r=d/2;%radius of the fracture,m
cden=0.05;%fracture density

ep=0.0001;
ht=[0.000001 0.01];
Ang=[0+ep,0+ep,pi/6,pi/4,pi/3,pi/2-ep];
MA=length(Ang);
f=-3:0.1:5;
f=10.^f;
Df=length(f);
v=zeros(Df,MA);
inQ=zeros(Df,MA);

for j=1:MA

ang=Ang(j);

if(j==1)
h=ht(1);
else
h=ht(2);
end

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;
den=dens*(1-por)+denf*por;

for FN=1:Df

FN

ft=f(FN);
w=2*pi*ft;
 
abk=(2*pi)/d*5;
DM=5000;
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
v(FN,j)=(2*pi*ft)/real(keff);
inQ(FN,j)=2*imag(keff)/real(keff);

end

end

%solutions of Galvin and Gurevich (2009) and Guo and Gurevich (2020b)

F=-3:0.5:5;
F=10.^F;
Df=length(F);
Dh=length(ht);
g=Gm/L;
b=vis/perm;
v1=zeros(Df,Dh);
inQ1=zeros(Df,Dh);
n0=cden/r^3;

for HN=1:Dh
    
    h=ht(HN);% when h tends to zero, Guo and Gurevich (2020b) reduces to Galvin and Gurevich (2009)

for FN=1:Df

FN
ft=F(FN);
w=2*pi*ft;
k0=10^(-3)*sqrt(den*w^2/L);%1/m
k1=10^(-3)*sqrt(den*w^2/H);%1/m
k2=sqrt(1i*w*b*H/(L*M));%1/m
k3=10^(-3)*sqrt(den*w^2/Gm);%1/m

abk=(2*pi)/d*5;
DM=5000;
y=abk/DM:abk/DM:abk;

%y=0.01:0.01:200;

YN=length(y);

q1=zeros(1,YN);
q2=zeros(1,YN);
q3=zeros(1,YN);

for lk=1:YN
    
yy=y(lk);
if yy<=k1
    q1(lk)=-1i*sqrt(k1^2-yy^2);
else
    q1(lk)=sqrt(yy^2-k1^2);
end
    q2(lk)=-1i*sqrt(k2^2-yy^2);
if yy<=k3
    q3(lk)=-1i*sqrt(k3^2-yy^2);
else
    q3(lk)=sqrt(yy^2-k3^2);
end

end

%T11=(k0^2-k1^2)*((M*L*h)./(q1*H)-2*Kf/k2^2*(1+(h*y.^2)./(2*q1))).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H;
T11=(k0^2-k1^2)*((M*L*h)./(q1*H)-2*Kf/k2^2).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H;
T12=(4*a*g*y.^2).*(y.^2-q2.*q3)-k2^2*(2*y.^2-k3^2);
%T21=q2.*(k2^2*(2*Kf/k2^2*(1+(h*y.^2)./(2*q2))-(M*L*h)./(q2*H)).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H);
T21=q2.*(k2^2*(2*Kf/k2^2-(M*L*h)./(q2*H)).*(2*y.^2-k3^2)+2*k3^2*Kf*a*M/H);
T22=2*a*g*(1-g)*k3^2*y;
T1=(T11.*T12)./(T21.*T22);        
T33=a*g*((2*y.^2-k3^2).^2-4*((y.^2).*q1).*q3)+(2*y.^2-k3^2)*(k0^2-k1^2);
T44=2*a*g*(1-g)*k3^2*q1.*y;
T2=T33./T44;
T=(1+(a*M*k3^2)./(H*(2*y.^2-k3^2))).*(T1-T2)-1;

N=50;
ds=r/N;
Ma=zeros(N,N);
V=zeros(N,1);
for i=1:N
    for j=1:N
       si=(i-1)*ds;
       sj=(j-1)*ds;
       
       TS=(T.*sin(si*y)).*sin(sj*y);
       
       My=trapz(y,TS)*2/pi;
       if(i==j)
           Ma(i,j)=My*ds+1;
       else
           Ma(i,j)=My*ds;
       end
       
    end
    
    p0=H-a*M;
    V(i)=-p0*si;
end

P=Ma\V;
sum=0;
for i=1:N
    si=(i-1)*ds;
    sum=sum+si*P(i)*ds;
end

sum=sum*2/pi;
g=Gm/L;
cof=-2*a*Gm*(1-g)/(L*(1-a*M/H));
A1=sum/cof;
f0=a*k1^2/L*A1;
ke=k1*(1+(2*pi*n0)/k1^2*f0);
v1(FN,HN)=w/real(ke);
inQ1(FN,HN)=2*imag(ke)/real(ke);

end

end

figure(1)
semilogx(f,v(:,1),'m','linewidth',2.5);
hold on;
semilogx(F,v1(:,1),'mo','MarkerSize',8);
hold on;
semilogx(f,v(:,2),'r','linewidth',2.5);
hold on;
semilogx(F,v1(:,2),'ro','MarkerSize',8);
hold on;
semilogx(f,v(:,3),'g','linewidth',2.5);
hold on;
semilogx(f,v(:,4),'b','linewidth',2.5);
hold on;
semilogx(f,v(:,5),'k','linewidth',2.5);
hold on;
semilogx(f,v(:,6),'c','linewidth',2.5);

xlabel('Frequency (Hz)','FontSize',15);
ylabel('{\itP}-wave velocity (m/s)','FontSize',15);
legend('Vp(0^{\circ}), {\it{c}}\rightarrow0','Galvin and Gurevich (2009)','Vp(0^{\circ})','Guo and Gurevich (2020b)','Vp(30^{\circ})','Vp(45^{\circ})','Vp(60^{\circ})','Vp(90^{\circ})');
set(gca, 'FontSize', 15);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

figure(2)
loglog(f,inQ(:,1),'m','linewidth',2.5);
hold on;
loglog(F,inQ1(:,1),'mo','MarkerSize',8);
hold on;
loglog(f,inQ(:,2),'r','linewidth',2.5);
hold on;
loglog(F,inQ1(:,2),'ro','MarkerSize',8);
hold on;
loglog(f,inQ(:,3),'g','linewidth',2.5);
hold on;
loglog(f,inQ(:,4),'b','linewidth',2.5);
hold on;
loglog(f,inQ(:,5),'k','linewidth',2.5);
hold on;
loglog(f,inQ(:,6),'c','linewidth',2.5);


xlabel('Frequency (Hz)','FontSize',15);
ylabel('Q_{p}^{-1}','FontSize',15);
legend('1/Qp(0^{\circ}), {\it{c}}\rightarrow 0','Galvin and Gurevich (2009)','1/Qp(0^{\circ})','Guo and Gurevich (2020b)','1/Qp(30^{\circ})','1/Qp(45^{\circ})','1/Qp(60^{\circ})','1/Qp(90^{\circ})');
set(gca, 'FontSize', 15);
xlim([0.001 100000]);
ylim([0.000001 0.1]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);








