function [f,g] = Ener_n_Ecc_withgrad( alp, dex, dey, dE, orb_param)
global mu
wEcc = 0.5;
wE = 0.5;
u = [-sin(alp) cos(alp)];

Ec = (orb_param(4)^2+orb_param(5)^2-1)*mu^2/2/orb_param(1)^2;

ex_end = orb_param(4) + sum(sum(dex.*u));
ey_end = orb_param(5) + sum(sum(dey.*u));
E_end  = Ec           + sum(sum(dE.*u));

f = wEcc*(ex_end^2+ey_end^2)+wE*(-0.5-E_end)^2;

if nargout > 1 % gradient required
    curl_u = [-cos(alp) -sin(alp)];
    a = sum(dex.*curl_u,2);
    b = sum(dey.*curl_u,2);
    c = sum(dE.*curl_u,2);
    g = wEcc*2*(ex_end)*a + wEcc*2*(ey_end)*b + wE*2*(0.5+E_end)*c;
end

end