function [f,g] = Ecc_n_H_3D(ctrl, dh, dhx, dhy, dex, dey, orb_param)
global wh_x we_x whxy_x

wh = 0.5;%wh_x;%0.66;
we = 0.1;%we_x;%0.04;
whxy = 0.4;%whxy_x;%0.35;
act_seg = length(dh(:,1));
alp = ctrl(1:act_seg);
bet = ctrl(act_seg+1:act_seg*2);

fr = -sin(alp).*cos(bet);
fn =  cos(alp).*cos(bet);
fh =  sin(bet);
frfh = fr.*fh;
fnfh = fn.*fh;
fhfh = fh.*fh;

u = [fr fn fh frfh fnfh fhfh];



h_end  = orb_param(1) + sum(sum(dh .*u));
hx_end = orb_param(2) + sum(sum(dhx.*u));
hy_end = orb_param(3) + sum(sum(dhy.*u));
ex_end = orb_param(4) + sum(sum(dex.*u));
ey_end = orb_param(5) + sum(sum(dey.*u));

f = wh*(1-h_end)^2+whxy*(hx_end^2+hy_end^2)+we*(ex_end^2+ey_end^2);

if nargout > 1 % gradient required
    curl_alp = [-cos(alp).*cos(bet)     -sin(alp).*cos(bet)      zeros(act_seg,1) ...
                -cos(alp).*cos(bet).*fh -sin(alp).*cos(bet).*fh  zeros(act_seg,1)];
    curl_bet = [ sin(alp).*sin(bet)     -cos(alp).*sin(bet)      cos(bet) ...
                -sin(alp).*(cos(bet).^2-sin(bet).^2) ...
                 cos(alp).*(cos(bet).^2-sin(bet).^2) ...
                 2*sin(bet).*cos(bet)];
    
    
    
    a1 = sum(dh.* curl_alp,2);
    b1 = sum(dhx.*curl_alp,2);
    c1 = sum(dhy.*curl_alp,2);
    d1 = sum(dex.*curl_alp,2);
    e1 = sum(dey.*curl_alp,2);
    
    a2 = sum(dh.* curl_bet,2);
    b2 = sum(dhx.*curl_bet,2);
    c2 = sum(dhy.*curl_bet,2);
    d2 = sum(dex.*curl_bet,2);
    e2 = sum(dey.*curl_bet,2);
    
    g_alp =   wh*2*(h_end-1)*a1 ...
           +whxy*2*(hx_end*b1 + hy_end*c1) ...
           +  we*2*(ex_end*d1 + ey_end*e1);
    g_bet =   wh*2*(h_end-1)*a2 ...
           +whxy*2*(hx_end*b2 + hy_end*c2) ...
           +  we*2*(ex_end*d2 + ey_end*e2);
    g = [g_alp; g_bet];

end

end
