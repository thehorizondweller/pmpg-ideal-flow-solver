function psi = calcPsi(x, y, a, theta)
    % Calculate p
    p = 1 - (x.*x)/(2*a)^2 - (y.*y)/(2*a)^2;

    % Calculate psi
    n = length(x);
    psi = zeros(1,n);
    for i = 1:n
        psi(i) = asinh(y(i)/(2*a*sin(theta(i))));
    end
end
