function a = conjug(c)
    % Discrete conjugation using complex FFT
    n = length(c);
    n1 = n / 2;
    k = 2:n1;
    a = fft(c);
    a(1) = 0;
    a(n1 + 1) = 0;
    a(k) = -1i * a(k);
    a(n1 + k) = 1i * a(n1 + k);
    a = ifft(a);
end
