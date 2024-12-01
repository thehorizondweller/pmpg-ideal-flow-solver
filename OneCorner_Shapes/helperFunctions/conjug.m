function a = conjug(c)
    % Discrete conjugation using complex FFT
    n = length(c);
    n1 = (length(c)) / 2;
    k = 2:n1;
    a = fft(c);
    a(1) = 0;
    a(n1 + 1) = 0;
    a(k) = -i * a(k);
    a(n1 + k) = i * a(n1 + k);
    a = ifft(a);
end
