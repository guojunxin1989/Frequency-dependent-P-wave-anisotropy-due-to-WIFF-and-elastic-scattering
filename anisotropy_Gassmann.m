function Cs=anisotropy_Gassmann(C,por,Ks,Kf)
% The function is used to calcuated the saturated stiffness matrix using
% anisotropic Gassmann (orthotropy). The parameters are as follows:
% C: dry stiffness matrix.
% por: porosity of the rock.
% Ks: bulk modulus of the grains.
% Kf: fluid modulus.

% The following are the three Biot's coefficients 
a1=1-(C(1,1)+C(1,2)+C(1,3))/(3*Ks);
a2=1-(C(1,2)+C(2,2)+C(2,3))/(3*Ks);
a3=1-(C(1,3)+C(2,3)+C(3,3))/(3*Ks);
Ma=(C(1,1)+C(2,2)+C(3,3)+2*C(1,2)+2*C(1,3)+2*C(2,3))/9;
M=Ks^2/(Ks*(1+por*(Ks/Kf-1))-Ma);%Biot's modulus 
%the followings are the effective stiffness of the saturated fractured rock
Cs11=C(1,1)+M*a1^2;
Cs12=C(1,2)+M*a1*a2;
Cs13=C(1,3)+M*a1*a3;
Cs22=C(2,2)+M*a2^2;
Cs23=C(2,3)+M*a2*a3;
Cs33=C(3,3)+M*a3^2;
Cs66=C(6,6);
Cs55=C(5,5);
Cs44=C(4,4);
Cs=[Cs11 Cs12 Cs13 0 0 0;
    Cs12 Cs22 Cs23 0 0 0;
    Cs13 Cs23 Cs33 0 0 0;
    0 0 0 Cs44 0 0;
    0 0 0 0 Cs55 0;
    0 0 0 0 0 Cs66];
end