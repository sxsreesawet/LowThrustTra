function flag = ...
    chkStop( h,hx,hy,ex,ey,tol_inc,tol_ecc,tol_a )
global mu

p = h^2/mu;

%==============================================
ecc = sqrt(ex^2+ey^2);
if ecc < (tol_ecc)
    flag_ecc = 1;
else
    flag_ecc = 0;
end
%==============================================
a = p/(1-ecc^2);
if (1-tol_a) < a && a < (1+tol_a)
    flag_a = 1;
else
    flag_a = 0;
end
%==============================================
i = asin(sqrt(hx^2+hy^2)/h)/pi*180;
if i < tol_inc
    flag_inc = 1;
else
    flag_inc = 0;
end
%==============================================


if flag_a && flag_ecc && flag_inc
    flag = 1;
else
    flag = 0;
end
end

