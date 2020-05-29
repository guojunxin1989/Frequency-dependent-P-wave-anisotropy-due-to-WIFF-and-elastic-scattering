clear;
close;
clc;

stype=1;%1:dry, 2:water saturated, 3: glycerin saturated

por=0.3460;%porosity of the background medium
dend=1.712;%dry background density
Ks=16.6;%bulk modulus of the solid
if(stype==1)
vp=2.529;
vs=1.558;
Gm=vs^2*dend;%shear modulus of the background medium
Lm=vp^2*dend;%P-wave modulus of the background medium,GPa
Km=Lm-4/3*Gm;
Kf=0.0001;%bulk modulus of the fluid
den=dend;%sample density
denf=0.0001;%fluid density
vis=0.00182*10^(-12);%viscosity of fluid, Gpa*s
elseif(stype==2)
Km=7.90;
Gm=3.942;
Kf=2.16;%bulk modulus of the fluid
denf=1;%fluid density
den=dend+denf*por;%sample density
vis=1*10^(-12);%viscosity of fluid, Gpa*s
end
perm=500*10^(-15);%permeability of background, m^2
d=5.5*10^(-3);%fracture diameter,m
r=d/2;% fracture radius, m
h=0.02*10^(-3);%thickness,m
cden=0.1;%fracture density

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;

f=[100000 1000];
tmp=0.0000001;
Ang=tmp:pi/50:pi/2+tmp;
MA=length(Ang);
v=zeros(MA,2);
inQ=zeros(MA,2);

for fn=1:2

w=2*pi*f(fn);
for j=1:MA
    
j

ang=Ang(j);
 
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
v(j,fn)=(2*pi*f(fn))/real(keff);
inQ(j,fn)=2*imag(keff)/real(keff);

end

end


N=2*MA-1;
Atp=tmp:pi/50:pi+tmp;
v1=zeros(N,2);
inQ1=zeros(N,2);

for k=1:2

for i=1:N
  
  if(i<=MA)
      j=i;
  else
      j=MA-(i-MA);
  end
  
  v1(i,k)=v(j,k);
  inQ1(i,k)=inQ(j,k);
end

end

if(stype==1)
%dry sample
A_lab=[10.9574 33.9603 56.461 78.9226 101.825 123.874 146.553 168.697];
v_lab=[1782.03 1914.47 2215.34 2467.1 2473.22 2205.7 1931.13 1782.9];
elseif(stype==2)
%water saturated
A_lab=[11.741 33.89 56.6072 78.6797 101.495 123.158 146.055 168.292];
v_lab=[2380.21 2423.29 2549.17 2688.79 2690.54 2547.37 2445.67 2378.41];
end
Re=0.006;
Rev=v_lab*Re;

Atp=Atp*180/pi;
  
figure(1)
plot(Atp,v1(:,1),'r','linewidth',1.5);
hold on;
plot(Atp,v1(:,2),'b','linewidth',1.5);
hold on;
plot(A_lab,v_lab,'ks','markersize',6,'MarkerFaceColor','k');

xlabel('Incidence angle','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('Predictions at 10^{5} Hz','Predictions at 10^{3} Hz','Experimental data')
set(gca, 'FontSize', 12);
xlim([0 180]);
set(gca,'XTick',[0 30 60 90 120 150 180]);

figure(2)
semilogy(Atp,inQ1(:,1),'r','linewidth',1.5);
hold on;
semilogy(Atp,inQ1(:,2),'b','linewidth',1.5);

xlabel('Incidence angle','FontSize',12);
ylabel('Q_{p}^{-1}','FontSize',12);
set(gca, 'FontSize', 12);
xlim([0 180]);
set(gca,'XTick',[0 30 60 90 120 150 180]);








