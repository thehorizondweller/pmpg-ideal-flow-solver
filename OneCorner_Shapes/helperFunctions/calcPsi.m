function psi = calcPsi(x, y, a, theta)
    % Calculate p
    p = 1 - (x.*x)/(2*a)^2 - (y.*y)/(2*a)^2;

    % Calculate psi
    n = length(x);
    psi = zeros(1,n);
    for i = 1:n
        % sqrt always spit out a positive value
        % Hence asinh will only spit positive values
        % psi(i) = asinh( sqrt( (p(i) + sqrt(p(i)^2 + y(i)^2/a^2)) / 2) );
        % if (x(i)^2 + y(i)^2) < a
        %     psi(i) = -psi(i);
        % end
        % psi(i) = y(i)/(2*a*sin(theta(i)));
        psi(i) = asinh(y(i)/(2*a*sin(theta(i))));
    end
end
