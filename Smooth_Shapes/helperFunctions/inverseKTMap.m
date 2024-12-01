function z = inverseKTMap(zeta, zeta1, zeta2, z1, z2, beta)
    % Inverse Transform from zeta to Z domain 
    q = ((zeta-zeta1)./(zeta-zeta2)).^beta;
    z = (z1-z2*q)./(1-q);
end
