function [ cond, result ] = main_optim( ini_val )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

time_start = cputime;
%%
%======================= Input and Setting ================================
% Initial Orbit
ini_R_perigee  = ini_val.orbit(1);          
ini_R_apogee   = ini_val.orbit(2);
inclination = ini_val.orbit(3);
% Spacecraft Parameters
Mass0 = ini_val.sp_param(1);
force1 = ini_val.sp_param(4);
I_sp1 = ini_val.sp_param(5);
force2 = ini_val.sp_param(2);
I_sp2 = ini_val.sp_param(3);
% Transfer Configuration
eclipse = ini_val.traf_fig(1);
swt_alt = ini_val.traf_fig(2);
% Terminal Condition
tol_inc = ini_val.terminal(1);
tol_ecc = ini_val.terminal(2);
tol_a   = ini_val.terminal(3);
max_rev = ini_val.terminal(4);
% Optimization Parameters
slices = ini_val.optim_param(1);               
Type = ini_val.optim_param(2);                     
option = ini_val.opt_alg;
%======================= End Input ========================================
%%
points = slices+1;
[h,hx,hy,ex,ey,m] = ...
    Initialize(Mass0,points,ini_R_perigee,ini_R_apogee,inclination);
% break;
global theta DU TU FU mu % SU MU
global Orbit_State Contr Isp Thr Trans_time
rev = 0;
% osc_count = 0;
while true
    rev = rev + 1
    Orbit_State(:,rev) = [h;hx;hy;ex;ey;m];
    Planar = chkPlanar(h,hx,hy,tol_inc);
    rp = h^2/mu/(1+sqrt(ex^2+ey^2));
    if Planar == 1
        if eclipse == 1
            ecl = chkEclipse(Orbit_State(1:5,rev));
        elseif eclipse == 0
            ecl = zeros(1,length(theta));
        elseif eclipse == 2
            ecl = chkEclipse(Orbit_State(1:5,rev))*2;
        else
            error('eclipse has to be 0 1 or 2');
        end
        
        if rp < (swt_alt+6378)/DU
            Isp = I_sp1/TU;
            Thr = force1/FU;
        else
            Isp = I_sp2/TU;
            Thr = force2/FU;
        end
        
        [h,hx,hy,ex,ey,m,alp,Detail(1:2,rev)] = ...
            Optim_planar(Thr,Isp,ecl,Orbit_State(:,rev),option,Type);
        
        Contr(:,:,rev) = [alp zeros(length(theta)-1,1)];
        Detail(3,rev) = 0.5*(1-h)^2+0.5*(ex^2+ey^2);
    elseif Planar ==0
        if eclipse == 1
            ecl = chkEclipse(Orbit_State(:,rev));
        elseif eclipse == 0
            ecl = zeros(1,length(theta));
        elseif eclipse == 2
            ecl = chkEclipse(Orbit_State(:,rev))*2;
        else
            error('eclipse has to be 0 1 or 2');
        end
        
        if rp <swt_alt/DU
            Isp = I_sp1/TU;
            Thr = force1/FU;
        else
            Isp = I_sp2/TU;
            Thr = force2/FU;
        end
%         error('Under Construction');
[h,hx,hy,ex,ey,m,ctrl,Detail(1:2,rev)] = ...
            Optim_3D(Thr,Isp,ecl,Orbit_State(:,rev),option,Type);
        Contr(:,:,rev) = ctrl;
        Detail(3,rev) = 0.333*(1-h)^2+0.333*(ex^2+ey^2)+0.333*(hx^2+hy^2);
    end
    
    
    stop = chkStop(h,hx,hy,ex,ey,tol_inc,tol_ecc,tol_a);
    if stop == 1
        display('Terminal Condition Satisfied! (^-'')9');
        Trans_time = sum(Orbit_State(1,:).^3/mu^2.*...
            sqrt(1./(1-Orbit_State(4,:).^2-Orbit_State(5,:).^2).^3));
        Orbit_State(:,rev+1) = [h;hx;hy;ex;ey;m];
        Final_mass = m*Mass0;
        break;
    elseif rev >= max_rev
        Orbit_State(:,rev+1) = [h;hx;hy;ex;ey;m];
        display('Beyond the maximum revolution! (-*-)');
        break;
    elseif rev>1
        if Detail(1,rev)>Detail(1,rev-1)
            display('Oscillation Detected! (-_-'')');
            h  = Orbit_State(1,rev);
            hx = Orbit_State(2,rev);
            hy = Orbit_State(3,rev);
            ex = Orbit_State(4,rev);
            ey = Orbit_State(5,rev);
            m  = Orbit_State(6,rev);
            rev = rev-1;
            
            if Type == 1
                Trans_time = sum(Orbit_State(1,:).^3/mu^2.*...
                sqrt(1./(1-Orbit_State(4,:).^2-Orbit_State(5,:).^2).^3));
                display('Fail to solve oscillation! (T^T)');
                break;
            else
                Type = 1;
            end
        end
    end
    
end
time_end = cputime;
computational_time = time_end-time_start

cond = 0;
result = 0;
end

