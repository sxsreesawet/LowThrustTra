function Rotation = Ro_rnh2dash(theta)
% =============== Error ===================================================
if ~isscalar(theta)
    error('theta has to be scalar')
end
% ================ Code ===================================================
Rotation = [cos(theta) -sin(theta) 0;
            sin(theta)  cos(theta) 0;
            0 0 1];

end

