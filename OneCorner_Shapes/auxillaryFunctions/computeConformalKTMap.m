function interDomData = computeConformalKTMap(z2,zetaPrime1,zetaPrime2,n,a,z,iterMax,accuracy,q)
    %% Step 2: Applying the KT Map
    % Get the angles at each corner
    angles = zeros(1,n);
    for i = 1:n
        prev = i-1;
        next = i+1;
        if i == 1
            prev = n;
        elseif i == n
            next = 1;
        end
        l1 = z(prev) - z(i);
        l2 = z(next) - z(i);
        angles(i) = acos((real(l1)*real(l2)+imag(l1)*imag(l2))/(abs(l1)*abs(l2)));
    end
    
    fig = figure('Visible','off');
    hold on;
    tempXAxis = linspace(1,n,n);
    plot(tempXAxis, angles*180/pi,"LineWidth",1.5);
    ylabel("Interior Angle (in Degrees");
    xlabel("Points on the Airfoil Surface (Starts on the TE and move CCW for increasingn point position)");
    title("Variation of Interior Angle with Point Location");
    hold off;
    outputFile = fullfile('./Results','02_Variation_of_Turning_Angle.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Get z1,z2,zeta1,zeta2 for the map
    [minAngle, minIdx] = min(angles);
    z1 = z(minIdx);

    % USER INPUTS - To be manually changed for debugging only
    % z2 = 0.01;
    % zetaPrime1 = 1;
    % zetaPrime2 = 0;

    beta = 2 - minAngle/pi;
    delta=1/beta;
    
    % Obtain the intermediate domain Zeta' by KT Map
    % Scales the zeta' to have z/zeta'->1 at infinity
    zetaPrime = directKTMap(z,z1,z2,zetaPrime1,zetaPrime2,delta);
    
    % Extract Real and Imag Parts
    xiPrime  = real(zetaPrime);
    etaPrime = imag(zetaPrime);
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    markerSize = 4; % Size of markers
    markerColor = 'b'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(xiPrime, etaPrime, markerSize, markerColor, 'filled', markerShape);
    xlabel("$$\xi'$$","Interpreter","latex");
    ylabel("$$\eta'$$","Interpreter","latex");
    title("$$\zeta'$$ Domain","Interpreter","latex");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','03_Intermediate_Domain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Centre the ZetaPrime domain
    zetaPrimeOffset = mean(zetaPrime);
    zetaPrimeCentered = zetaPrime - zetaPrimeOffset;

    %% Step 3: Calculate Theta and Psi

    theta = zeros(1,n);
    for i=1:n
        theta(i) = angle(zetaPrimeCentered(i));
        if theta(i) < 0
            theta(i) = 2*pi + theta(i);
        end
    end
    disp("Step 2: Theta for Intermediate Domain computed.");
    
    psi = log(abs(zetaPrimeCentered)/a);
    disp("Step 3: Psi for Intermediate Domain computed.");
    
    
    % Variation of Psi with Theta
    thetaPsi = [theta', psi']; % Zip the values of theta and psi
    sorted_thetaPsi = sortrows(thetaPsi,1); % Sort the theta values
    theta = sorted_thetaPsi(:,1)';
    psi = sorted_thetaPsi(:,2)';
    
    %% Step 4: Interpolation
    
    %Make Periodic Extensions of Theta and Psi
    thetaExt = [theta-4*pi theta-2*pi theta theta+2*pi theta+4*pi];
    psiExt = [psi psi psi psi psi];
    
    % Interpolation Method
    theta_interp = linspace(0,2*pi,n*100);
    psi_interp = interp1(thetaExt, psiExt, theta_interp, 'spline'); % 'spline' interpolation method
    
    fig = figure('Visible','off');
    subplot(2,1,1);
    hold on;
    grid on;
    plot(theta, psi, "LineWidth",2,"Color","b");
    xlabel("$$\theta$$","Interpreter","latex");
    ylabel("$$\psi$$","Interpreter","latex");
    title("$$\psi(\theta)$$","Interpreter","latex");
    % axis equal;
    hold off;
    subplot(2,1,2);
    hold on;
    grid on;
    plot(theta_interp, psi_interp, "LineWidth",2,"Color","r");
    xlabel("$$\theta$$","Interpreter","latex");
    ylabel("$$\psi$$ (Interpolated)","Interpreter","latex");
    title("Periodically Interpolated $$\psi(\theta)$$","Interpreter","latex");
    % axis equal;
    hold off;
    outputFile = fullfile('./Results','04_Radial_vs_Angular_Distribution.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % WARNING !!! - n gets updated here due to interpolation
    n = length(theta_interp);
    
    disp("Step 4: Interpolation of Psi and Theta complete.");
    
    %% Step 5 - Calculating Epsilon 
    f = a*exp(psi_interp).*exp(1i*theta_interp);

    % USER INPUTS - To be manually changed for debugging only
    % iterMax = 1e04;
    % accuracy = 1e-08;

    theta_phi = theoext(n, f, iterMax, accuracy, psiExt, thetaExt, a);
    
    phi = 2*pi*(0:n-1)/n;
    epsilon = phi-theta_phi;
    fig = figure('Visible','off');
    hold on;
    grid on;
    plot(phi, theta_phi,"LineWidth",2,'Color','k');
    xlabel("$$\phi$$","Interpreter","latex");
    ylabel("$$\theta$$","Interpreter","latex");
    title("Boundary Correspondence Function");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','05_Boundary_Correspondence_Function.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    plot(phi, epsilon,'LineWidth',2,'Color','r');
    xlabel("$$\phi$$","Interpreter","latex");
    ylabel("$$\varepsilon$$","Interpreter","latex");
    title("Angular Distortion");
    % axis equal;
    hold off;
    outputFile = fullfile('./Results','06_Angular_Distortion_Function.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    
    %% STEP 5B - Calculate Psi as a function of Phi
    psi_phi = interp1(thetaExt,psiExt,theta_phi,'spline');
    % Plot Psi as a function of Phi
    fig = figure('Visible','off');
    hold on;
    grid on;
    plot(phi, psi_phi, "LineWidth",2,"Color","k");
    xlabel("$$\phi$$","Interpreter","latex");
    legend({"$$\psi(\phi)$$"},"Interpreter","latex");
    title("$$\psi(\phi)$$","Interpreter","latex");
    hold off;
    outputFile = fullfile('./Results','07_Psi_vs_Phi.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 5B: Psi as a function of Phi obtained.");
    
    
    %% STEP 6 - Calculate coefficients of the Zeta->Z' Transformation
    psi_0 = mean(psi_phi);
    r = a*exp(psi_0);

    % USER INPUTS - Manually update for debugging only
    % q = 100; % Order upto which the series expansion is to be considered for now

    A = zeros(1,q); % Real Part of the Coefficients
    B = zeros(1,q); % Imag Part of the Coefficients
    for i = 1:q
        % Finding the fourier coefficients
        A(i) = ((r^i)/pi) * trapz(phi, psi_phi.*cos(i*phi));
        B(i) = ((r^i)/pi) * trapz(phi, psi_phi.*sin(i*phi));
    
    end
    
    % Retrace Zeta' and Z starting from Zeta to validate
    zeta = a * exp(psi_0) * exp(1i*phi); % Points on the circle in the zeta domain
    zetaPrimeRetrace = zeros(1,n);
    for i = 1:n
        expSum = 0;
        for m = 1:q
            expSum = expSum + (A(m)+1i*B(m))/zeta(i)^m;
        end
        zetaPrimeRetrace(i) = zeta(i)*exp(expSum);
    end
    
    % Reversing the Centering Operation
    zetaPrimeRetrace = zetaPrimeRetrace + zetaPrimeOffset;
    % Parameters are the same as declared before for directKTMap
    zRetrace = inverseKTMap(zetaPrimeRetrace, zetaPrime1, zetaPrime2, z1, z2, beta);
    
    fig = figure('Visible','off');
    
    subplot(1,3,1)
    hold on;
    grid on;
    markerSize = 4; % Size of markers
    markerColor = 'b'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zeta),imag(zeta),markerSize,markerColor,'filled',markerShape);
    xlabel("$$\xi$$","Interpreter","latex");
    ylabel("$$\eta$$","Interpreter","latex");
    title("$$\zeta Domain$$", "Interpreter","latex");
    axis equal;
    hold off;
    
    subplot(1,3,2)
    hold on;
    grid on;
    markerSize = 10; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zetaPrimeRetrace),imag(zetaPrimeRetrace),markerSize,markerColor,'filled',markerShape);
    hold on;
    markerSize = 3;
    markerColor = 'r';
    markerShape = 's';
    scatter(real(zetaPrime),imag(zetaPrime),markerSize,markerColor,"filled",markerShape);
    xlabel("$$\xi'$$","Interpreter","latex");
    ylabel("$$\eta'$$","Interpreter","latex");
    legend(["$$\zeta'$$ (Retraced)","$$\zeta'$$ (Original)"],"Interpreter","latex");
    title("$$\zeta' Domain$$", "Interpreter","latex");
    axis equal;
    hold off;
    
    subplot(1,3,3)
    hold on;
    grid on;
    markerSize = 10; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zRetrace),imag(zRetrace),markerSize,markerColor,'filled',markerShape);
    hold on;
    markerSize = 2; % Size of markers
    markerColor = 'r'; % Color of markers (e.g., 'r' for red)
    markerShape = 's'; % Shape of markers (e.g., 's' for square)
    scatter(real(z),imag(z),markerSize,markerColor,'filled',markerShape);
    xlabel("$$x$$","Interpreter","latex");
    ylabel("$$y$$","Interpreter","latex");
    legend(["$$z$$ (Retraced)","$$z$$ (Original)"],"Interpreter","latex");
    title("$$z Domain$$", "Interpreter","latex");
    axis equal;
    hold off;
    
    jointTitle = sprintf("Order of Fourier Coefficients used: %d",q);
    sgtitle(jointTitle);

    outputFile = fullfile('./Results','08_All_Three_Domains.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 6: Mapping Complete!");

    interDomData.minAngle = minAngle;
    interDomData.z1 = z1;
    interDomData.beta = beta;
    interDomData.A = A;
    interDomData.B = B;
    interDomData.psi_0 = psi_0;
    interDomData.zetaPrimeOffset = zetaPrimeOffset;
    interDomData.n = n;
    interDomData.zRetrace = zRetrace;
    interDomData.zeta = zeta;
end