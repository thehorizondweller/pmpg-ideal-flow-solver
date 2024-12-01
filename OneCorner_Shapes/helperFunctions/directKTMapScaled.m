function zeta = directKTMapScaled(z, z1, z2, zeta1, zeta2, delta)
    % Karman–Trefftz from z to zeta plane
    w = ((z(1)-z1)/(z(1)-z2))^delta;
    zeta(1) = (zeta1 - zeta2*w)./(1-w);
    w = ((z(2)-z1)/(z(2)-z2))^delta;
    zeta(2) = (zeta1 - zeta2*w)./(1-w);
    for j=3:length(z)
        branch = w;
        w = ((z(j)-z1)/(z(j)-z2));
        w = crootbr(w,branch,delta);
        zeta(j) = (zeta1 - zeta2*w)/(1-w);
    end
    zeta = zeta*(delta*((z1-z2)/(zeta1-zeta2)));
end