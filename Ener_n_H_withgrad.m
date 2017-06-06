function [f,g] = Ener_n_H_withgrad( alp, dh, dE, orb_param)
global mu
wh = 0.5;
wE = 0.5;
u = [-sin(alp) cos(alp)];

Ec = (orb_param(4)^2+orb_param(5)^2-1)*mu^2/2/orb_param(1)^2;

h_end  = orb_param(1) + sum(sum(dh.*u));
E_end  = Ec           + sum(sum(dE.*u));

f = wh*(1-h_end)^2+wE*(-0.5-E_end)^2;

if nargout > 1 % gradient required
    curl_u = [-cos(alp) -sin(alp)];
    a = sum(dh.*curl_u,2);
    b = sum(dE.*curl_u,2);
    g = wh*2*(h_end-1)*a + wE*2*(0.5+E_end)*b;
end

end
