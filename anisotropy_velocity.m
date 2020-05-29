function [Vp,Vsv,Vsh]=anisotropy_velocity(C11,C12,C22,C44,C55,C66,ang,den)
% This function is used to calculate the velocities in the 1-2 plane. 
% The velocities are dependent on the incident angles. The parameters are
% as follows:
% C11, C12, C22, and C66 are the stiffness coefficients in the 1-2 plane.
% ang: The angle here is the angle relative to the 1-axis.
% den: density of the rock. 
% This function can also be used to calculate the velocities in other planes, 
% just need to replace the above coefficients with the coefficients of the
% corresponding plane.

Ms=(C66*sin(ang)^2+C11*cos(ang)^2)*(C22*sin(ang)^2+C66*cos(ang)^2)-(C12+C66)^2*sin(ang)^2*cos(ang)^2;
Cof=sqrt((C66+C22*sin(ang)^2+C11*cos(ang)^2)^2-4*Ms);
Vp=sqrt((C66+C22*sin(ang)^2+C11*cos(ang)^2+Cof)/(2*den));
Vsv=sqrt((C66+C22*sin(ang)^2+C11*cos(ang)^2-Cof)/(2*den));
Vsh=sqrt((C55*cos(ang)^2+C44*sin(ang)^2)/den);
end