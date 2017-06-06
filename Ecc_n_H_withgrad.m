function [f,g] = Ecc_n_H_withgrad( alp, dh, dex ,dey, orb_param)
global wh_x we_x
wh = 0.75;%0.75;%wh_x;%0.5;
we = 0.25;%0.25;%we_x;%0.5;
u = [-sin(alp) cos(alp)];

h_end  = orb_param(1) + sum(sum(dh.*u));
ex_end = orb_param(4) + sum(sum(dex.*u));
ey_end = orb_param(5) + sum(sum(dey.*u));

f = wh*(1-h_end)^2+we*(ex_end^2+ey_end^2);

if nargout > 1 % gradient required
    curl_u = [-cos(alp) -sin(alp)];
    a = sum(dh.*curl_u,2);
    b = sum(dex.*curl_u,2);
    c = sum(dey.*curl_u,2);
    g = wh*2*(h_end-1)*a + we*2*(ex_end)*b + we*2*(ey_end)*c;
end

end

