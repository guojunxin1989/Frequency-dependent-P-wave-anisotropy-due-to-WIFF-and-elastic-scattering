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
h=0.01;%thickness,m
cden=0.05;%fracture density

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;
den=dens*(1-por)+denf*por;

f=[0.001,0.1,300,3000];
Ang=0:pi/100:pi;
Df=length(f);
MA=length(Ang);
v=zeros(Df,MA);
inQ=zeros(Df,MA);

for j=1:MA

ang=Ang(j);

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

% calculate P-wave disperion and attenuation using Guo et al.ги2018a,2018b)

Kmf=0.00000001;
Gmf=0.00000001;
g=Gm/L;
a1=1;
a2=1;
h=0.01;
a3=h/a1;
pf=4/3*pi*a3*cden;
opt=2;

Ang1=0:pi/20:pi;
N=length(Ang1);
Cl=zeros(1,6);
Ch=zeros(1,6);
v1=zeros(3,N);
inQ1=zeros(3,N);

%low frequency limit
Sf=General_Eshelby_model(Km,Gm,Kmf,Gmf,a1,a2,a3,pf,opt);
ZN=Sf(3,3);
ZT=Sf(4,4);
La=L-2*Gm;
C0=[L,La,La,0,0,0;
    La,L,La,0,0,0;
    La,La,L,0,0,0;
    0,0,0,Gm,0,0;
    0,0,0,0,Gm,0;
    0,0,0,0,0,Gm;
    ];
S0=inv(C0);
S1=S0+Sf;
C1=inv(S1);
por1=por*(1-pf)+pf;
Cs1=anisotropy_Gassmann(C1,por1,Ks,Kf);
Cl(1)=Cs1(3,3);
Cl(2)=Cs1(2,3);
Cl(3)=Cs1(2,2);
Cl(4)=Cs1(6,6);
Cl(5)=Cs1(5,5);
Cl(6)=Cs1(4,4);
den=dens*(1-por)+denf*por;

%intermediate frequencies
Kms=Km+a^2*M;
Sf=General_Eshelby_model(Kms,Gm,Kmf,Gmf,a1,a2,a3,pf,opt);
La=H-2*Gm;
C0=[H,La,La,0,0,0;
    La,H,La,0,0,0;
    La,La,H,0,0,0;
    0,0,0,Gm,0,0;
    0,0,0,0,Gm,0;
    0,0,0,0,0,Gm;
    ];
S0=inv(C0);
S1=S0+Sf;
C1=inv(S1);
Cs2=anisotropy_Gassmann(C1,pf,Kms,Kf);
Ch(1)=Cs2(3,3);
Ch(2)=Cs2(2,3);
Ch(3)=Cs2(2,2);
Ch(4)=Cs2(6,6);
Ch(5)=Cs2(5,5);
Ch(6)=Cs2(4,4);

Cm=(Ch-Cl)./Cl;

T=(2*(H-a*M)^2*(2-4*a*g+3*a^2*g^2)*r^2*cden*vis)/(15*Gm*g*(1-g)^2*H*L*perm);
Sc=(pi*cden)/r;
C2=Cs2(3,3);
C1=Cs1(3,3);
Lc=pf/ZN;
Gc=pf/ZT;
Kc=Lc-4/3*Gc;
ac=1-Kc/Ks;
porf=1;
Mc=((ac-porf)/Ks+porf/Kf)^(-1);
Hc=Lc+ac^2*Mc;
G=(2*Sc*C2*(a*M/H-ac*Mc/Hc)^2)/sqrt((M*L*vis)/(H*perm));
F1=((C2-C1)/(C1*G))^2;
F2=(C2-C1)^3/(2*C2*C1^2*T*G^2);

f=[0.001 0.1 300];
Df=length(f);

for i=1:Df
 
ft=f(i);
s=(1./Ch).*(1+Cm/(1-F2+F2*sqrt(1-1i*(2*pi*ft*F1)/F2^2)));
c=1./s;

for j=1:N
ang=Ang1(j);
[Vp,Vsv,Vsh]=anisotropy_velocity(c(1),c(2),c(3),c(4),c(5),c(6),ang,den);
v1(i,j)=(real(1/Vp))^(-1)*1000;
inQ1(i,j)=-imag(Vp^2)/real(Vp^2);
end

end

Ang1=Ang1*180/pi;
Ang=Ang*180/pi;

figure(1)
plot(Ang,v(1,:),'r','linewidth',1.5);
hold on;
plot(Ang1,v1(1,:),'ro','MarkerSize',6);
hold on;
plot(Ang,v(2,:),'g','linewidth',1.5);
hold on;
plot(Ang1,v1(2,:),'go','MarkerSize',6);
hold on;
plot(Ang,v(3,:),'b','linewidth',1.5);
hold on;
plot(Ang1,v1(3,:),'bo','MarkerSize',6);
hold on;
plot(Ang,v(4,:),'k','linewidth',1.5);

xlabel('Incidence angle','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('{\itf} = 0.001 Hz','LF validation','{\itf} = 0.1 Hz','CF validation','{\itf} = 300 Hz','IF validation','{\itf} = 3000 Hz');
set(gca, 'FontSize', 12);
xlim([0 180]);
set(gca,'XTick',[0 30 60 90 120 150 180]);

figure(2)
semilogy(Ang,inQ(1,:),'r','linewidth',1.5);
hold on;
semilogy(Ang1,inQ1(1,:),'ro','MarkerSize',6);
hold on;
semilogy(Ang,inQ(2,:),'g','linewidth',1.5);
hold on;
semilogy(Ang1,inQ1(2,:),'go','MarkerSize',6);
hold on;
semilogy(Ang,inQ(3,:),'b','linewidth',1.5);
hold on;
semilogy(Ang1,inQ1(3,:),'bo','MarkerSize',6);
hold on;
semilogy(Ang,inQ(4,:),'k','linewidth',1.5);

xlabel('Incidence angle','FontSize',12);
ylabel('Q_{p}^{-1}','FontSize',12);
legend('{\itf} = 0.001 Hz','LF validation','{\itf} = 0.1 Hz','CF validation','{\itf} = 300 Hz','IF validation','{\itf} = 3000 Hz');
set(gca, 'FontSize', 12);
xlim([0 180]);
ylim([0.000001 0.1]);
%ylim([0 0.05]);
set(gca,'XTick',[0 30 60 90 120 150 180]);







