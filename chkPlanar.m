function planar = chkPlanar( h,hx,hy,tol_inc )
i = asin(sqrt(hx^2+hy^2)/h)/pi*180;
if i < tol_inc
    planar = 1;
else
    planar = 0;
end
end

