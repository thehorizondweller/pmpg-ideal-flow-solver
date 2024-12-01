function theta = calcTheta(x, y, a)
    % Calculate p
    p = 1 - (x.*x)/(4*a*a) - (y.*y)/(4*a*a);

    % Calculate theta
    n = length(x);
    theta = zeros(1,n);
    for i = 1:n
        soln = asin( sqrt( (p(i) + sqrt(p(i)^2 + y(i)^2/a^2)) / 2) );
        if soln == 0
            if x(i) > 0
                theta(i) = 0;
            else
                theta(i) = pi;
            end
        else
            if x(i) > 0 && y(i) > 0
                theta(i) = soln;
            elseif x(i) < 0 && y(i) > 0 
                theta(i) = pi - soln;
            elseif x(i) < 0 && y(i) < 0
                theta(i) = pi + soln;
            elseif x(i) > 0 && y(i) < 0
                theta(i) = 2*pi - soln;
            end
        end
    end
