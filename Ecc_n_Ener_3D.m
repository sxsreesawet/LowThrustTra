function [f,g] = Ecc_n_Ener_3D(ctrl, dE, dhx, dhy, dex, dey, orb_param)
global mu

wE = 0.3;
we = 0.08;
whxy = 0.62;
act_seg = length(dE(:,1));
alp = ctrl(1:act_seg);
bet = ctrl(act_seg+1:act_seg*2);

fr = -sin(alp).*cos(bet);
fn =  cos(alp).*cos(bet);
fh =  sin(bet);
frfh = fr.*fh;
fnfh = fn.*fh;
fhfh = fh.*fh;

u = [fr fn fh frfh fnfh fhfh];

Ec = (orb_param(4)^2+orb_param(5)^2-1)*mu^2/2/orb_param(1)^2;

E_end  = Ec + sum(sum(dE .*u));
hx_end = orb_param(2) + sum(sum(dhx.*u));
hy_end = orb_param(3) + sum(sum(dhy.*u));
ex_end = orb_param(4) + sum(sum(dex.*u));
ey_end = orb_param(5) + sum(sum(dey.*u));

f = wE*(-0.5-E_end)^2+whxy*(hx_end^2+hy_end^2)+we*(ex_end^2+ey_end^2);

if nargout > 1 % gradient required
    curl_alp = [-cos(alp).*cos(bet)     -sin(alp).*cos(bet)      zeros(act_seg,1) ...
                -cos(alp).*cos(bet).*fh -sin(alp).*cos(bet).*fh  zeros(act_seg,1)];
    curl_bet = [ sin(alp).*sin(bet)     -cos(alp).*sin(bet)      cos(bet) ...
                -sin(alp).*(cos(bet).^2-sin(bet).^2) ...
                 cos(alp).*(cos(bet).^2-sin(bet).^2) ...
                 2*sin(bet).*cos(bet)];
    
    
    
    a1 = sum(dE.* curl_alp,2);
    b1 = sum(dhx.*curl_alp,2);
    c1 = sum(dhy.*curl_alp,2);
    d1 = sum(dex.*curl_alp,2);
    e1 = sum(dey.*curl_alp,2);
    
    a2 = sum(dE.* curl_bet,2);
    b2 = sum(dhx.*curl_bet,2);
    c2 = sum(dhy.*curl_bet,2);
    d2 = sum(dex.*curl_bet,2);
    e2 = sum(dey.*curl_bet,2);
    
    g_alp =   wE*2*(0.5+E_end)*a1 ...
           +whxy*2*(hx_end*b1 + hy_end*c1) ...
           +  we*2*(ex_end*d1 + ey_end*e1);
    g_bet =   wE*2*(0.5+E_end)*a2 ...
           +whxy*2*(hx_end*b2 + hy_end*c2) ...
           +  we*2*(ex_end*d2 + ey_end*e2);
    g = [g_alp; g_bet];

end

end
