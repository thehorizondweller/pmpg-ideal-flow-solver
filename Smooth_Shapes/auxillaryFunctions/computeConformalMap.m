function conformalMapData = computeConformalMap(a,psi,theta,psiExt,thetaExt,n,z,zetaPrime,q)
    psi_interp = psi;
    theta_interp = theta;
    %% Step 5 - Calculating Epsilon 
    f = a*exp(psi_interp).*exp(1i*theta_interp);
    iterMax = 2e03;
    accuracy = 1e-08;
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
    outputFile = fullfile('./Results','04_Boundary_Correspondence_Function.png');
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
    outputFile = fullfile('./Results','05_Angular_Distortion.png');
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
    outputFile = fullfile('./Results','06_Psi_vs_Phi.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 5B: Psi as a function of Phi obtained.");
    
    
    %% STEP 6 - Calculate coefficients of the Zeta->Z' Transformation
    psi_0 = mean(psi_phi);
    r = a*exp(psi_0);
    % q = 100; - User input provided; Change manually for debugging only
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
    
    zRetrace2 = zetaPrimeRetrace + a^2./zetaPrimeRetrace;
    
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

    markerSize = 10;
    markerColor = 'r';
    markerShape = 's';
    scatter(real(zetaPrime),imag(zetaPrime),markerSize,markerColor,"filled",markerShape);
    hold on;
    markerSize = 2; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zetaPrimeRetrace),imag(zetaPrimeRetrace),markerSize,markerColor,'filled',markerShape);
    
    xlabel("$$\xi'$$","Interpreter","latex");
    ylabel("$$\eta'$$","Interpreter","latex");
    legend(["$$\zeta'$$ (Original)","$$\zeta'$$ (Retraced)"],"Interpreter","latex");
    title("$$\zeta' Domain$$", "Interpreter","latex");
    axis equal;
    hold off;
    
    subplot(1,3,3)
    hold on;
    grid on;
    markerSize = 10; % Size of markers
    markerColor = 'r'; % Color of markers (e.g., 'r' for red)
    markerShape = 's'; % Shape of markers (e.g., 's' for square)
    scatter(real(z),imag(z),markerSize,markerColor,'filled',markerShape);
    hold on;
    markerSize = 2; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zRetrace2),imag(zRetrace2),markerSize,markerColor,'filled',markerShape);
   
    xlabel("$$x$$","Interpreter","latex");
    ylabel("$$y$$","Interpreter","latex");
    legend(["$$z$$ (Original)","$$z$$ (Retraced)"],"Interpreter","latex");
    title("$$z Domain$$", "Interpreter","latex");
    axis equal;
    hold off;
    
    jointTitle = sprintf("Order of Fourier Coefficients used: %d",q);
    sgtitle(jointTitle);

    outputFile = fullfile('./Results','07_ConformalMap_Three_Domains.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 6: Mapping Complete!");

    conformalMapData.A = A;
    conformalMapData.B = B;
    conformalMapData.psi_0 = psi_0;
end