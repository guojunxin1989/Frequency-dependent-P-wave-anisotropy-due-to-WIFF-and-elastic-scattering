function H=General_Eshelby_model(Km,Gm,Kmf,Gmf,a1,a2,a3,por,opt)
% This function is used to calculate the excess compliance matrix of
% fractures which is filled by any elastic materials in the isotropic
% elastic background matrix, Mura, 1987, section 11
% Km, Gm: the bulk and shear moduli of the background matrix 
% Kmf,Gmf: the bulk and shear moduli of the crack inclusion 
% a1, a2, a3: the three radius of the ellipsoidal crack along the x1,x2,and x3 axis
% por: crack porosity
% opt: shape of the inclusions:
% 1: sphere; 2: penny shaped cracks a1=a2>>a3; 3: oblate spheroid a1=a2>a3;
% 4: elliptic cylinder a3=infinity 5:Flat ellipsoid a1>a2>>a3; 6:prolate
% spheroid a1>a2=a3; 7: general case a1>a2>a3
% H: the obtained fracture compliance matrix
% Copyright: Junxin Guo

C0=[Km+4/3*Gm,Km-2/3*Gm,Km-2/3*Gm,0,0,0;
    Km-2/3*Gm,Km+4/3*Gm,Km-2/3*Gm,0,0,0;
    Km-2/3*Gm,Km-2/3*Gm,Km+4/3*Gm,0,0,0;
    0,0,0,Gm,0,0;
    0,0,0,0,Gm,0;
    0,0,0,0,0,Gm;
    ];
C1=C0;
for i=1:6
    for j=4:6
        C1(i,j)=2*C0(i,j);
    end
end

J=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0;
   0,0,0,0.5,0,0;
   0,0,0,0,0.5,0;
   0,0,0,0,0,0.5;
   ];

v=(3*Km-2*Gm)/(2*(3*Km+Gm));

if(opt==1)
%   S11=(7-5*v)/(15*(1-v));
%   S12=(5*v-1)/(15*(1-v));
%   S44=(4-5*v)/(15*(1-v));
a2=a1;
a3=a1;
I1=4*pi/3;
I2=I1;
I3=I1;
I11=4*pi/(5*a1^2);
I22=I11;
I33=I11;
I12=I11;
I13=I11;
I23=I11;
elseif(opt==2)
%   S11=(13-8*v)/(32*(1-v))*pi*(a3/a1);
%   S22=S11;
%   S33=1-(1-2*v)/(1-v)*pi/4*(a3/a1);
%   S12=(8*v-1)/(32*(1-v))*pi*(a3/a1);
%   S21=S12;
%   S13=(2*v-1)/(8*(1-v))*pi*(a3/a1);
%   S23=S13;
%   S31=v/(1-v)*(1-(4*v+1)/(8*v)*pi*(a3/a1));
%   S32=S31;
%   S66=(7-8*v)/(32*(1-v))*pi*(a3/a1);
%   S55=1/2*(1+(v-2)/(1-v)*pi/4*(a3/a1));
%   S44=S55;
I1=pi^2*a3/a1;
I2=I1;
I3=4*pi-2*pi^2*a3/a1;
I12=3*pi^2*a3/(4*a1^3);
I13=3*(4/3*pi-pi^2*a3/a1)/a1^2;
I23=I13;
I11=3*pi^2*a3/(4*a1^3);
I22=I11;
I33=4/3*pi/a3^2;
elseif(opt==3)
  I1=2*pi*a1^2*a3/(a1^2-a3^2)^(3/2)*(acos(a3/a1)-a3/a1*(1-a3^2/a1^2)^(1/2));
  I2=I1;
  I3=4*pi-2*I1;
  I13=(I1-I3)/(a3^2-a1^2);
  I23=I13;
  I33=(4*pi/a3^2-2*I13)/3;
  I12=pi/a1^2-1/4*I13;
  I11=I12;
  I22=I12;
  
elseif(opt==4)
    
   I1=4*pi*a2/(a1+a2);
   I2=4*pi*a1/(a1+a2);
   I3=0;
   
   I12=4*pi/(a1+a2)^2;
   I11=(4*pi/a1^2-I12)/3;
   I22=(4*pi/a2^2-I12)/3;
   I13=I1/a3^2;
   I23=I2/a3^2;
   I33=0;
elseif(opt==5)
   k=(a1^2-a2^2)/a1^2;
   e=@(x)(1-k*sin(x).^2).^(1/2);
   E=integral(@(x)e(x),0,pi/2);
   f=@(x)(1-k*sin(x).^2).^(-1/2);
   F=integral(@(x)f(x),0,pi/2);
   I1=4*pi*a2*a3*(F-E)/(a1^2-a2^2);
   I2=4*pi*a3*E/a2-4*pi*a2*a3*(F-E)/(a1^2-a2^2);
   I3=4*pi-4*pi*a3*E/a2;
   I12=(4*pi*a3*E/a2-8*pi*a2*a3*(F-E)/(a1^2-a2^2))/(a1^2-a2^2);
   I23=(4*pi-8*pi*a3*E/a2+4*pi*a2*a3*(F-E)/(a1^2-a2^2))/a2^2;
   I13=(4*pi-4*pi*a2*a3*(F-E)/(a1^2-a2^2)-4*pi*a3*E/a2)/a1^2;
   I33=4*pi/(3*a3^2);
   I11=(4*pi/a1^2-I12-I13)/3;
   I22=(4*pi/a2^2-I12-I23)/3;
elseif(opt==6)
   I2=(2*pi*a1*a3^2)/(a1^2-a3^2)^(3/2)*((a1/a3)*(a1^2/a3^2-1)^(1/2)-acosh(a1/a3));
   I3=I2;
   I1=4*pi-2*I2;
   I12=(I2-I1)/(a1^2-a2^2);
   I11=(4*pi/a1^2-2*I12)/3;
   I22=(4*pi/a2^2-(I2-I1)/(a1^2-a2^2))/4;
   I33=I22;
   I23=I22;
   I13=(4*pi/a1^2-3*I11-I12);
elseif(opt==7)
   ang=asin((1-a3^2/a1^2)^(1/2));
   k=((a1^2-a2^2)/(a1^2-a3^2))^(1/2);
   e=@(x)(1-k^2*sin(x).^2).^(1/2);
   E=integral(@(x)e(x),0,ang);
   f=@(x)(1-k^2*sin(x).^2).^(-1/2);
   F=integral(@(x)f(x),0,ang);
   I1=(4*pi*a1*a2*a3)/((a1^2-a2^2)*(a1^2-a3^2)^(1/2))*(F-E);
   I3=(4*pi*a1*a2*a3)/((a2^2-a3^2)*(a1^2-a3^2)^(1/2))*(a2*(a1^2-a3^2)^(1/2)/(a1*a3)-E);
   I2=4*pi-I1-I3;
   I12=(I2-I1)/(a1^2-a2^2);
   I13=(I3-I1)/(a1^2-a3^2);
   I23=(I3-I2)/(a2^2-a3^2);
   I11=(4*pi/a1^2-I12-I13)/3;
   I22=(4*pi/a2^2-I12-I23)/3;
   I33=(4*pi/a3^2-I13-I23)/3;
end

I21=I12;
I31=I13;
I32=I23;

S11=3/(8*pi*(1-v))*a1^2*I11+(1-2*v)/(8*pi*(1-v))*I1;
S22=3/(8*pi*(1-v))*a2^2*I22+(1-2*v)/(8*pi*(1-v))*I2;
S33=3/(8*pi*(1-v))*a3^2*I33+(1-2*v)/(8*pi*(1-v))*I3;
  
S12=1/(8*pi*(1-v))*a2^2*I12-(1-2*v)/(8*pi*(1-v))*I1;
S21=1/(8*pi*(1-v))*a1^2*I21-(1-2*v)/(8*pi*(1-v))*I2;
S13=1/(8*pi*(1-v))*a3^2*I13-(1-2*v)/(8*pi*(1-v))*I1;
S31=1/(8*pi*(1-v))*a1^2*I31-(1-2*v)/(8*pi*(1-v))*I3;
S23=1/(8*pi*(1-v))*a3^2*I23-(1-2*v)/(8*pi*(1-v))*I2;
S32=1/(8*pi*(1-v))*a2^2*I32-(1-2*v)/(8*pi*(1-v))*I3;
  
S44=(a2^2+a3^2)/(16*pi*(1-v))*I23+(1-2*v)/(16*pi*(1-v))*(I2+I3);
S55=(a1^2+a3^2)/(16*pi*(1-v))*I13+(1-2*v)/(16*pi*(1-v))*(I1+I3);
S66=(a1^2+a2^2)/(16*pi*(1-v))*I12+(1-2*v)/(16*pi*(1-v))*(I1+I2);
  
S=[S11,S12,S13,0,0,0;
   S21,S22,S23,0,0,0;
   S31,S32,S33,0,0,0;
   0,0,0,S44,0,0;
   0,0,0,0,S55,0;
   0,0,0,0,0,S66;
   ];

Q=C1*(J-S);

M0=inv(C0);

C2=[Kmf+4/3*Gmf,Kmf-2/3*Gmf,Kmf-2/3*Gmf,0,0,0;
    Kmf-2/3*Gmf,Kmf+4/3*Gmf,Kmf-2/3*Gmf,0,0,0;
    Kmf-2/3*Gmf,Kmf-2/3*Gmf,Kmf+4/3*Gmf,0,0,0;
    0,0,0,Gmf,0,0;
    0,0,0,0,Gmf,0;
    0,0,0,0,0,Gmf;
    ];
M1=inv(C2);

M2=M0;
M2(4,4)=1/4*M0(4,4);
M2(5,5)=1/4*M0(5,5);
M2(6,6)=1/4*M0(6,6);

M3=M1;
M3(4,4)=1/4*M1(4,4);
M3(5,5)=1/4*M1(5,5);
M3(6,6)=1/4*M1(6,6);

DM=M3-M2;
C3=DM\J;

C4=C3;
for j=1:6
    for i=4:6
        C4(i,j)=C3(i,j)/2;
    end
end

C=C4+Q;
H1=C\J;

H=H1;
for j=1:6
    for i=4:6
        H(i,j)=H1(i,j)/2;
    end
end

H=por*H;

for i=1:6
    for j=1:6
        if(i<4&&j<4)
            H(i,j)=H(i,j);
        elseif(i>3&&j>3)
            H(i,j)=4*H(i,j);
        else
            H(i,j)=2*H(i,j);
        end
    end
end

end


   


