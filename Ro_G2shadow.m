function Rot = Ro_G2shadow( Day )
%   Return a rotation matrix from earth fixed to shadow fixed frame

% Earth tilted angle of rotation
tilt = 28.5/180*pi;
R1 = [1 0 0;
    0 cos(tilt) sin(tilt);
    0 -sin(tilt) cos(tilt)];
% Assuming Northern Hemisphere equinox happens on 0.00 am March 22.
Do = 79;
delay = 79+180;
ang = 180 + ((Day-Do+delay)/365.25*360);
ang_rad = ang/180*pi;
R3 = [cos(ang_rad) sin(ang_rad) 0;
    -sin(ang_rad) cos(ang_rad) 0;
    0 0 1];
Rot = R3*R1;

end

