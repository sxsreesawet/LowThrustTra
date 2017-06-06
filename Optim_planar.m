function [h_f,hx_f,hy_f,ex_f,ey_f,mass,alpha,detail] = ...
                                Optim_planar(F,Isp,ecl,state,options,type )
global mu theta g0
h = state(1);
ex = state(4);
ey = state(5);
m = state(6);

a = ex*sin(theta.') - ey*cos(theta.');
b = 1+ex*cos(theta.')+ey*sin(theta.');
gh = h^5/mu^3/m*F*[zeros(length(theta),1) 1./b.^3];
gex = h^4/mu^2/m*F*[sin(theta.')./b.^2  2*cos(theta.')./b.^2+sin(theta.').*a./b.^3];
gey = h^4/mu^2/m*F*[-cos(theta.')./b.^2 2*sin(theta.')./b.^2-cos(theta.').*a./b.^3];
gE = h^2/mu/m*F*[a./(b.^2) 1./b];

gMass = h^3/mu^2./b.^2;
act_seg = 0;
for i = 1:(length(theta)-1) % trapezoidal rule
     if ecl(i) == 0 && ecl(i+1) == 0
            act_seg = act_seg+1;
            dh(act_seg,:)  = (theta(i+1)-theta(i))*(gh(i,:)+gh(i+1,:))/2;
            dex(act_seg,:) = (theta(i+1)-theta(i))*(gex(i,:)+gex(i+1,:))/2;
            dey(act_seg,:) = (theta(i+1)-theta(i))*(gey(i,:)+gey(i+1,:))/2;
            dmass(act_seg) = (theta(i+1)-theta(i))*(gMass(i) + gMass(i+1))/2 ...
                             *F/Isp/g0;
            dE(act_seg,:)  = (theta(i+1)-theta(i))*(gE(i,:)+gE(i+1,:))/2;
            
     elseif ecl(i) == 2 || ecl(i+1) == 2
            act_seg = act_seg+1;
            dh(act_seg,:)  = (theta(i+1)-theta(i))*(gh(i,:)+gh(i+1,:))/2/4;
            dex(act_seg,:) = (theta(i+1)-theta(i))*(gex(i,:)+gex(i+1,:))/2/4;
            dey(act_seg,:) = (theta(i+1)-theta(i))*(gey(i,:)+gey(i+1,:))/2/4;
            dmass(act_seg) = (theta(i+1)-theta(i))*(gMass(i) + gMass(i+1))/2/4 ...
                             *F/Isp/g0;
            dE(act_seg,:)  = (theta(i+1)-theta(i))*(gE(i,:)+gE(i+1,:))/2/4;
     end
end
alp0 = zeros(act_seg,1);

if type == 1
    [alp_res,fval,exitflag,output] = ...
            fminunc(@(alp)Ecc_n_H_withgrad(alp,dh,dex,dey,state)...
            ,alp0,options);
elseif type == 2
    [alp_res,fval,exitflag,output] = ...
            fminunc(@(alp)Ener_n_H_withgrad(alp,dh,dE,state)...
            ,alp0,options);
elseif type == 3
    [alp_res,fval,exitflag,output] = ...
            fminunc(@(alp)Ener_n_Ecc_withgrad(alp,dex,dey,dE,state)...
            ,alp0,options);
else
    Error('Type has to be 1 or 2');
end

f_star = [-sin(alp_res) cos(alp_res)];

h_f = h + sum(sum(dh.*f_star));
ex_f = ex + sum(sum(dex.*f_star));
ey_f = ey + sum(sum(dey.*f_star));

hx_f = h_f/h*state(2);
hy_f = h_f/h*state(3);
mass = state(6)-sum(dmass);

count = 0;
for i = 1:(length(theta)-1)
    if ecl(i) == 0 && ecl(i+1) == 0
        count = count+1;
        alpha(i,1) = alp_res(count);
    elseif ecl(i) == 2 || ecl(i+1) == 2
        count = count+1;
        alpha(i,1) = alp_res(count);
    else
        alpha(i,1) = inf;
    end
end
detail = [fval;exitflag];
end

