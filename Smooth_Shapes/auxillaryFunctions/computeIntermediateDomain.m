function interDomData = computeIntermediateDomain(x,y,z,a,n)
    %% Step 2: Calculate Theta    
    theta = calcTheta(x, y, a);
    disp("Step 2: Theta for Non-offset Intermediate Domain computed.");

    %% Step 3: Calculate Psi
    psi = calcPsi(x, y, a, theta);
    nan_ind_psi = find(isnan(psi));
    psi(nan_ind_psi) = psi(nan_ind_psi+1);
    while ~isempty(nan_ind_psi)
        nan_ind_psi = find(isnan(psi));
        psi(nan_ind_psi) = psi(nan_ind_psi+1);
    end
    disp("Step 3: Psi for Non-offset Intermediate Domain computed.");
    
    % Variation of Psi with Theta
    thetaPsi = [theta', psi']; % Zip the values of theta and psi
    sorted_thetaPsi = sortrows(thetaPsi,1); % Sort the theta values
    theta = sorted_thetaPsi(:,1)';
    psi = sorted_thetaPsi(:,2)';
    
    %% Step 3B: Interpolation
    %Make Periodic Extensions of Theta and Psi
    thetaExt = [theta-4*pi theta-2*pi theta theta+2*pi theta+4*pi];
    psiExt = [psi psi psi psi psi];
    
    % Interpolation Method
    theta_interp = linspace(0,2*pi,n*100);
    psi_interp = interp1(thetaExt, psiExt, theta_interp, 'spline'); % 'spline' interpolation method
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    plot(theta, psi, "LineWidth",4,"Color","b");
    hold on;
    plot(theta_interp, psi_interp, "LineWidth",2,"Color","r");
    xlabel("$$\theta$$","Interpreter","latex");
    ylabel("$$\psi$$","Interpreter","latex");
    title("$$\psi(\theta)$$","Interpreter","latex");
    legend("Original Psi","Interpolated Psi")
    hold off;
    outputFile = fullfile('./Results','02_Radial_vs_Angular_Distribution.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % WARNING !!! - n gets updated here due to interpolation
    n = length(theta_interp);
    
    disp("Step 3B: Interpolation of Psi and Theta complete.");
    
    %% Step 4 - Z' plane
    zetaPrime = zeros(1,n);
    for i=1:n
        zetaPrime(i) = a * exp(psi_interp(i)) * exp(1i*theta_interp(i));
    end
    
    zRetrace1 = zetaPrime + a^2./zetaPrime;
    fig = figure('Visible','off');
    hold on;
    grid on;
    markerSize = 4; % Size of markers
    markerColor = 'b'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zetaPrime),imag(zetaPrime),markerSize,markerColor,'filled',markerShape);
    hold on;
    markerSize = 4; % Size of markers
    markerColor = 'r'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(zRetrace1),imag(zRetrace1),markerSize,markerColor,'filled',markerShape);
    hold on;
    markerSize = 4; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 's'; % Shape of markers (e.g., 's' for square)
    scatter(real(z),imag(z),markerSize,markerColor,'filled',markerShape);
    xlabel("X'");
    ylabel("Y'");
    legend("Z Prime Plane", "Z (Airfoil) Retraced", "Z (Airfoil) Original");
    title("Intermediate Z' and Z (Retraced and Original) Domains");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','03_IntermediateDomain.png');
    exportgraphics(fig, outputFile, 'Resolution',600);
    
    disp("Step 4: Intermediate Domain Computed.")

    interDomData.theta = theta_interp; % Angular distribution
    interDomData.psi = psi_interp; % Radial distribution
    interDomData.n = n; % Number of points on the contour after interpolation
    interDomData.zetaPrime = zetaPrime; % contour of the intermediate domain

    interDomData.thetaExt = thetaExt;
    interDomData.psiExt = psiExt;
end