clear;
close;
clc;

stype=2;%1:dry, 2:water saturated

por=0.33;%porosity of the background medium
Ks=37;%bulk modulus of the solid
if(stype==1)
Km=5.7188;%bulk modulus of the background medium,GPa
Gm=6.6758;%shear modulus of the background medium
Kf=0.0001;%bulk modulus of the fluid
den=1.736;%sample density
denf=0.0001;%fluid density
vis=0.00182*10^(-12);%viscosity of fluid, Gpa*s
elseif(stype==2)
Km=8.40;%bulk modulus of the background medium,GPa
Gm=7.45;%shear modulus of the background medium
Kf=2.19;%bulk modulus of the fluid
den=2.065;%sample density
denf=1;%fluid density
vis=1*10^(-12);%viscosity of fluid, Gpa*s
end
perm=21*10^(-15);%permeability of background, m^2
r=2.91*10^(-3);%radius of the fracture,m
d=r*2;%fracture diameter
as=0.0429;%fracture aspect ratio
h=d*as;%thickness,m
cden=0.0314;%fracture density

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;

f=1.5*100000;
w=2*pi*f;
tmp=0.0000001;
Ang=tmp:pi/50:pi/2+tmp;
MA=length(Ang);
v=zeros(MA,1);
inQ=zeros(MA,1);

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
v(j)=(2*pi*f)/real(keff);
inQ(j)=2*imag(keff)/real(keff);

end


N=2*MA-1;
Atp=tmp:pi/50:pi+tmp;
Atp=Atp*180/pi;
v1=zeros(N,1);
inQ1=zeros(N,1);

for i=1:N
  
  if(i<=MA)
      j=i;
  else
      j=MA-(i-MA);
  end
  
  v1(i)=v(j);
  inQ1(i)=inQ(j);
end

A_lab=[0 45 90 135 180];
if(stype==1)
%dry sample 
v_lab=[2631 2800 2902 2781 2637];
att_lab=[0.116148 0.018677 0 0.057782 0.116148];
elseif(stype==2)
%water saturated
v_lab=[3125 3182 3262 3202 3127];
att_lab=[0.0633574 0.0341602 0 0.0410284 0.0637394];
%att_lab=[0.127099 0.0979018 0.0637416 0.10477 0.127481];
else
end
Re=0.006;
Rev=v_lab*Re;
Rea=att_lab*0.2;

  
figure(1)

plot(Atp,v1,'r','linewidth',1.5);
hold on;
errorbar(A_lab,v_lab,Rev,'rs','markersize',6,'MarkerFaceColor','r');
hold on;

xlabel('Incidence angle','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('Theoretical prediction','experimental data');
set(gca, 'FontSize', 12);
xlim([0 180]);
set(gca,'XTick',[0 30 60 90 120 150 180]);

figure(2)

plot(Atp,inQ1,'r','linewidth',1.5);
hold on;
errorbar(A_lab,att_lab,Rea,'rs','markersize',6,'MarkerFaceColor','r');
hold on;


xlabel('Incidence angle','FontSize',12);
ylabel('Q_{p}^{-1}','FontSize',12);
legend('Theoretical prediction','experimental data');
set(gca, 'FontSize', 12);
xlim([0 180]);
set(gca,'XTick',[0 30 60 90 120 150 180]);










