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
perm0=10^(-15);%permeability of background, m^2
d=1;%diameter of the fracture,m
r=d/2;%radius of the fracture,m
h=0.01;%thickness,m
cden=0.05;%fracture density

Per=[perm0,10*perm0,100*perm0,1000*perm0];
DP=length(Per);
f=-3:0.1:5;
f=10.^f;
Df=length(f);
tpE=zeros(Df,DP);
tpD=zeros(Df,DP);

for PN=1:DP

perm=Per(PN);
L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;
den=dens*(1-por)+denf*por;

for FN=1:Df

FN

VA=zeros(1,3);
ep=0.001;
Ang=[0+ep,pi/4,pi/2-ep];

for VN=1:3

ang=Ang(VN);

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
VA(VN)=(2*pi*ft)/keff;
end

tpE(FN,PN)=(VA(3)-VA(1))/VA(1);
tpD(FN,PN)=4*(VA(2)/VA(1)-1)-tpE(FN,PN);

end

end

figure(1)
semilogx(f,real(tpE(:,1)),'r','linewidth',1.5);
hold on;
semilogx(f,real(tpE(:,2)),'g','linewidth',1.5);
hold on;
semilogx(f,real(tpE(:,3)),'b','linewidth',1.5);
hold on;
semilogx(f,real(tpE(:,4)),'c','linewidth',1.5);
hold on;
semilogx(f,real(tpD(:,1)),'r--','linewidth',1.5);
hold on;
semilogx(f,real(tpD(:,2)),'g--','linewidth',1.5);
hold on;
semilogx(f,real(tpD(:,3)),'b--','linewidth',1.5);
hold on;
semilogx(f,real(tpD(:,4)),'c--','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Velocity anisotropy parameters','FontSize',12);
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

figure(2)
semilogx(f,imag(tpE(:,1)),'r','linewidth',1.5);
hold on;
semilogx(f,imag(tpE(:,2)),'g','linewidth',1.5);
hold on;
semilogx(f,imag(tpE(:,3)),'b','linewidth',1.5);
hold on;
semilogx(f,imag(tpE(:,4)),'c','linewidth',1.5);
hold on;
semilogx(f,imag(tpD(:,1)),'r--','linewidth',1.5);
hold on;
semilogx(f,imag(tpD(:,2)),'g--','linewidth',1.5);
hold on;
semilogx(f,imag(tpD(:,3)),'b--','linewidth',1.5);
hold on;
semilogx(f,imag(tpD(:,4)),'c--','linewidth',1.5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Attenuation anisotropy parameters','FontSize',12);
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);








