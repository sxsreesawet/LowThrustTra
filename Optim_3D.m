function [h_f,hx_f,hy_f,ex_f,ey_f,mass,ctrl,detail] = ...
                                Optim_3D(F,Isp,ecl,state,options,type )
global mu theta g0
h = state(1);
hx = state(2);
hy = state(3);
ex = state(4);
ey = state(5);
m = state(6);


a = ex*sin(theta.')   -  ey*cos(theta.');
b = 1+ex*cos(theta.') + ey*sin(theta.');
C = Ro_dash2G([hx;hy;sqrt(h^2-hx^2-hy^2)]);
tanx = -hy/sqrt(h^2-hy^2);

dex1dt = -sin(theta.')*(-h/m);                                           %Correct!
dex2dt = 2*h/m*cos(theta.')+h./m.*a./b.*sin(theta.');                    %Correct!
dex3dt = tanx*h*-ey/mu/m*sin(theta.')./b;
dey1dt = cos(theta.')*(-h/m);                                            %Correct!
dey2dt = 2*h/m*sin(theta.')-h./m.*a./b.*cos(theta.');                    %Correct!
dey3dt = tanx*h*ex/mu/m*sin(theta.')./b;

dhdt  = h^2/mu/m*F*[zeros(length(theta),1) 1./b zeros(length(theta),1)]; %Correct!
dhxdt = h^2/mu/m*F* ...
        [zeros(length(theta),1), ...
        C(1,3)./b, ...
        C(1,1)*sin(theta.')./b-C(1,2)*cos(theta.')./b];
dhydt = h^2/mu/m*F* ...
        [zeros(length(theta),1), ...
        C(2,3)./b, ...
        C(2,1)*sin(theta.')./b-C(2,2)*cos(theta.')./b];
% dhxdt = h/m/mu*F*...
%         [zeros(length(theta),1), ...
%          hx./b, ...
%          h*sqrt(h^2-hx^2-hy^2)/sqrt(h^2-hy^2).*sin(theta.')./b+...
%          hx*hy/sqrt(h^2-hy^2).*cos(theta.')./b];
% dhydt = h/m/mu*F*...
%         [zeros(length(theta.'),1), ...
%         hy./b,...
%         -sqrt(h^2-hy^2).*cos(theta.')./b];
    
dexdt = [dex1dt dex2dt dex3dt]*F;
deydt = [dey1dt dey2dt dey3dt]*F;
dEnerdt = mu/h*F*[a b zeros(length(theta),1)];

dtdtheta = [h^3/mu^2./b.^2, ...
           -h^3/mu^2./b.^2.*sin(theta.')*tanx*h^4/mu^3./b.^3/m*F];
expand = [dtdtheta(:,1) dtdtheta(:,1) dtdtheta(:,1), ...
          dtdtheta(:,2) dtdtheta(:,2) dtdtheta(:,2)];
M_h  = [dhdt.*expand(:,1:3)  dhdt.*expand(:,4:6)];
M_hx = [dhxdt.*expand(:,1:3) dhxdt.*expand(:,4:6)];
M_hy = [dhydt.*expand(:,1:3) dhydt.*expand(:,4:6)];
M_ex = [dexdt.*expand(:,1:3) dexdt.*expand(:,4:6)];
M_ey = [deydt.*expand(:,1:3) deydt.*expand(:,4:6)];
M_Ener = [dEnerdt.*expand(:,1:3) dEnerdt.*expand(:,4:6)];

M_m =  dtdtheta;


act_seg = 0;
for i = 1:(length(theta)-1) % trapezoidal rule
     if ecl(i) == 0 && ecl(i+1) == 0
            act_seg = act_seg+1;
            dh(act_seg,:)  = (theta(i+1)-theta(i))*(M_h(i,:)+M_h(i+1,:))/2;
            dhx(act_seg,:) = (theta(i+1)-theta(i))*(M_hx(i,:)+M_hx(i+1,:))/2;
            dhy(act_seg,:) = (theta(i+1)-theta(i))*(M_hy(i,:)+M_hy(i+1,:))/2;
            dex(act_seg,:) = (theta(i+1)-theta(i))*(M_ex(i,:)+M_ex(i+1,:))/2;
            dey(act_seg,:) = (theta(i+1)-theta(i))*(M_ey(i,:)+M_ey(i+1,:))/2;
            dmass(act_seg,:) = (theta(i+1)-theta(i))*(M_m(i,:)+M_m(i+1,:))/2 ...
                             *F/Isp/g0;
            dE(act_seg,:)  = (theta(i+1)-theta(i))*(M_Ener(i,:)+M_Ener(i+1,:))/2;
     
     elseif ecl(i) == 2 || ecl(i+1) == 2 %ecl == 2 for 1/4 thrust
            act_seg = act_seg+1;
            dh(act_seg,:)  = (theta(i+1)-theta(i))*(M_h(i,:)+M_h(i+1,:))/2;
            dhx(act_seg,:) = (theta(i+1)-theta(i))*(M_hx(i,:)+M_hx(i+1,:))/2;
            dhy(act_seg,:) = (theta(i+1)-theta(i))*(M_hy(i,:)+M_hy(i+1,:))/2;
            dex(act_seg,:) = (theta(i+1)-theta(i))*(M_ex(i,:)+M_ex(i+1,:))/2;
            dey(act_seg,:) = (theta(i+1)-theta(i))*(M_ey(i,:)+M_ey(i+1,:))/2;
            dmass(act_seg,:) = (theta(i+1)-theta(i))*(M_m(i,:)+M_m(i+1,:))/2 ...
                             *F/Isp/g0;
            dE(act_seg,:)  = (theta(i+1)-theta(i))*(M_Ener(i,:)+M_Ener(i+1,:))/2;
            
            %Thrust Loss Adjustment
            dh(act_seg,:) = [dh(act_seg,1:3)/4 dh(act_seg,4:6)/16];
            dhx(act_seg,:) = [dhx(act_seg,1:3)/4 dhx(act_seg,4:6)/16];
            dhy(act_seg,:) = [dhy(act_seg,1:3)/4 dhy(act_seg,4:6)/16];
            dex(act_seg,:) = [dex(act_seg,1:3)/4 dex(act_seg,4:6)/16];
            dey(act_seg,:) = [dey(act_seg,1:3)/4 dey(act_seg,4:6)/16];
            dmass(act_seg,:) = [dmass(act_seg,1)/4 dmass(act_seg,2)/16];
            dE(act_seg,:) = [dE(act_seg,1:3)/4 dE(act_seg,4:6)/16];
            
     end
end

ctrl0 = zeros(act_seg*2,1);
if type == 1
    [ctrl_res,fval,exitflag,output] = ...
            fminunc(@(ctrl)Ecc_n_H_3D(ctrl,dh,dhx,dhy,dex,dey,state)...
            ,ctrl0,options);
elseif type == 2
    error('Type 2 for 3D trajectory is not available \_(-.-)_/');
elseif type == 3
%     error('Type 3 for 3D trajectory is not available \_(-.-)_/');
    [ctrl_res,fval,exitflag,output] = ...
            fminunc(@(ctrl)Ecc_n_Ener_3D(ctrl,dE,dhx,dhy,dex,dey,state)...
            ,ctrl0,options);
else
    error('Type has to be 1');
end

alp_star = ctrl_res(1:act_seg);
bet_star = ctrl_res(act_seg+1:act_seg*2);

fr_star = -sin(alp_star).*cos(bet_star);
fn_star =  cos(alp_star).*cos(bet_star);
fh_star =  sin(bet_star);
frfh_star = fr_star.*fh_star;
fnfh_star = fn_star.*fh_star;
fhfh_star = fh_star.*fh_star;

f_star = [fr_star fn_star fh_star frfh_star fnfh_star fhfh_star];

h_f  = h  + sum(sum(dh .*f_star));
hx_f = hx + sum(sum(dhx.*f_star));
hy_f = hy + sum(sum(dhy.*f_star));
ex_f = ex + sum(sum(dex.*f_star));
ey_f = ey + sum(sum(dey.*f_star));
mass = m  - sum(sum([dmass(:,1) dmass(:,2).*fh_star]));

count = 0;
for i = 1:(length(theta)-1)
    if ecl(i) == 0 && ecl(i+1) == 0
        count = count+1;
        ctrl(i,1:2) = [ctrl_res(count) ctrl_res(count+act_seg)];
    elseif ecl(i) == 2 || ecl(i+1) == 2
        count = count+1;
        ctrl(i,1:2) = [ctrl_res(count) ctrl_res(count+act_seg)];
    else
        ctrl(i,1:2) = [inf inf];
    end
end
detail = [fval;exitflag];
end

