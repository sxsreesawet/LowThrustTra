function Rotation = Ro_dash2G( H )
% =============== Error ===================================================
[m,n] = size(H);
if (m~=3)||(n~=1)
    error('Input has to be 3x1 vector ([hx;hy;hz])');
end
% ================ Code ===================================================
hx = H(1);
hy = H(2);
hz = H(3);

x_dash_hat = [hz;0;-hx]/sqrt(hx^2+hz^2);
h_hat = [hx;hy;hz]/sqrt(hx^2+hy^2+hz^2);
xn_dash_hat = [0         -h_hat(3) h_hat(2);
               h_hat(3)     0       -h_hat(1);
              -h_hat(2)   h_hat(1)    0]*x_dash_hat;
Rotation = [x_dash_hat xn_dash_hat h_hat];

end

