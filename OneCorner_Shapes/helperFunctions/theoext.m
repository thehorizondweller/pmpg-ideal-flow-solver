function theta_phi = theoext(divs, f, iterMax, accuracy, psiExt, thetaExt, a)
    phi = 2*pi*(0:divs-1)/divs;
    theta = phi;
    theta1 = theta;
    disp("Iteration No. || Error");
    error = 1e08;
    currIter = 1;
    while (error>accuracy) && (currIter<iterMax)
        c = log(abs(f));
        c = -conjug(c);
        theta = real(c) + phi;
        error = max(abs(theta-theta1));
        theta1 = theta;
        fprintf('%6.0f %30.14e \n', currIter, error);
        psiEval = interp1(thetaExt, psiExt, theta, 'spline');
        f = a*exp(psiEval).*exp(1i*theta);
        currIter = currIter + 1;
    end
    fprintf('Iterations to reach desired accuracy: %d \n', currIter);

    % for iter = 1:iterMax
    %     c = log(abs(f));
    %     c = -conjug(c);
    %     theta = real(c) + phi;
    %     error = max(abs(theta-theta1));
    %     theta1 = theta;
    %     fprintf('%6.0f %30.14e \n', iter, error);
    %     psiEval = interp1(thetaExt, psiExt, theta, 'spline');
    %     f = a*exp(psiEval).*exp(1i*theta);
    % end
    theta_phi = theta;
end

