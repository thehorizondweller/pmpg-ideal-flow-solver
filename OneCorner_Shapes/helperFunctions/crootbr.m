function zb = crootbr(z, branch, delta)
    % choosing continuous complex branch of delta root
    b1 = real(branch);
    b2 = imag(branch);
    z = z^delta;
    q = real(z);
    r = imag(z);
    if (q*b1+r*b2 < 0)
        zb = z * exp(-1i*2*pi*delta);
    else
        zb = z;
    end
end
