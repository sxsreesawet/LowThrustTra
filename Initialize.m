function [h0,hx0,hy0,ex0,ey0,m0] = Initialize(mass0,segments,rp,ra,i)
global DU TU SU MU FU mu

DU = 42164;                 %distance unit, Km
TU = sqrt(42164^3/398600);  %time unit, s
SU = DU/TU;                 %speed unit, Km/sec
MU = mass0;                    %mass Unit, Kg
FU = MU*DU*1000/TU^2;       %force unit(N)
mu = 1;

global rEarth g0 theta
rEarth = 6378/DU;           %normalized radius of the Earth
g0 = 9.81/1000/DU*TU^2;     %normalized gravitational acc
theta = linspace(0,2*pi,segments);

if rp ~= 0 || ra ~=0
    a = (rp+ra)/2/DU;
    ex0 = ra/DU/a-1;
    ey0 = 0;

    h0 = mu*sqrt(a*(1-ex0^2));
    hx0 = -sin(i/180*pi)*h0;
    hy0 = 0;
    m0 = 1;

else
    % J.T. Betts Initial Condition. 
    % Isp = 1849.34774 s, Mass0 = 1000 kg , Force = 1.445 N.
%     r = [-6878.14000000000;7.40253327326573e-13;-4.01924763248280e-13]/DU;
%     v = [-9.32275413600020e-16;-6.69008883015199;3.63242186141846]/DU*TU;

    % Robert D. Falck I.C. (1)
    % Isp = 3300 s, Mass0 = 1200 kg, Force = [0.401706386930699] N.
%     r = [6927;0;0]/DU;
%     v = [0;6.66645959243842;3.61959147841841]/DU*TU;
    
    % Robert D. Falck I.C. (2)
    % Isp = 1800 s, Mass0 = 1200 kg, Force = [0.311579953965478] N.
    r = [11359.0700000000;0;0]/DU;
    v = [0;9.00933471253942;4.89166861353448]/DU*TU;
    %======================================================================
    rv = [r;v];
    Param = RV2Param(rv)
    
    h0 = norm(Param(1:3));
    hx0 = Param(1);
    hy0 = Param(2);
    ex0 = Param(4);
    ey0 = Param(5);

    i0 = acos(Param(3)/h0)/pi*180
    m0 = 1;
    g0 = 9.80665000000000/1000/DU*TU^2;
end

