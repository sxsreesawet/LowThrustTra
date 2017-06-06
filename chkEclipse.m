function isEcl = chkEclipse( Orbit_State )
global theta rEarth mu
h = Orbit_State(1);
hx = Orbit_State(2);
hy = Orbit_State(3);
ex = Orbit_State(4);
ey = Orbit_State(5);

hz = sqrt(h^2-hx^2-hy^2);


r = h^2/mu./(1+ex*cos(theta)+ey*sin(theta));
r = [r; zeros(2,length(r))];

for i=1:length(theta)
    Rxyz = Ro_dash2G([hx;hy;hz])*Ro_rnh2dash(theta(i))*r(:,i);
    if Rxyz(1)>0 && sqrt(Rxyz(2)^2 + Rxyz(3)^2)<rEarth
        isEcl(1,i) = 1;
    else
        isEcl(1,i) = 0;
        
    end
end


end

