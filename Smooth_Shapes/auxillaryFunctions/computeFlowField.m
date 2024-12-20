function flowFieldData = computeFlowField(radialPts,angularPts,terminalR,a,psi_0,U,alpha,beta,kuttaNormGDev,A,B)
    
    %% STEP 7 - Build the Domains
    
    % Build the Domain

    % Converted to USER INPUTS - Manually update for debugging only
    % radialPts = 41;
    % angularPts = 101;
    % terminalR = 10;

    r = a*exp(psi_0); % Radius of the circle in the zeta domain
    domainGrids = flowDomainBuild(radialPts, angularPts, terminalR, a, r, A, B);
    ZETA = domainGrids.ZETA;
    ZETA_PRIME = domainGrids.ZETA_PRIME;
    Z = domainGrids.Z;
    
    % Plot the Domains
    fig = figure('Visible','off');
    subplot(2,2,1);
    hold on; grid on;
    markerSize = 3; % Size of markers
    markerColor = 'b'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(ZETA),imag(ZETA),markerSize,markerColor,'filled',markerShape);
    xlabel("$$\xi$$","Interpreter","latex");
    ylabel("$$\eta$$","Interpreter","latex");
    title("$$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(2,2,2);
    hold on; grid on;
    markerSize = 3; % Size of markers
    markerColor = 'r'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(ZETA_PRIME),imag(ZETA_PRIME),markerSize,markerColor,'filled',markerShape);
    xlabel("$$\xi '$$","Interpreter","latex");
    ylabel("$$\eta '$$","Interpreter","latex");
    title("$$\zeta '$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(2,2,3);
    hold on; grid on;
    markerSize = 3; % Size of markers
    markerColor = 'k'; % Color of markers (e.g., 'r' for red)
    markerShape = 'c'; % Shape of markers (e.g., 's' for square)
    scatter(real(Z),imag(Z),markerSize,markerColor,'filled',markerShape);
    xlabel("$$x$$","Interpreter","latex");
    ylabel("$$y$$","Interpreter","latex");
    title("$$z$$ Domain","Interpreter","latex");
    axis equal;
 
    hold off;
    outputFile = fullfile('./Results','08_Grid_Domains.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 7: Grid Constructed for all 3 domains.")
    
    %% STEP 8A - Velocity Field Computations in the Zeta Domain
    
    % Converted to USER INPUTS - Manually update for debugging only
    % U = 1; % Far Field Conditions
    % alpha = 5; % In degrees

    alpha = alpha * pi / 180; % Convert to Radians
    
    % Initial Guess for Gamma based on the Kutta Condition
    kutta_normG = sin(alpha+beta);
    G = 4*pi*U*r*kutta_normG;
    
    % Complex Potential in the Zeta Domain
    FZeta = U*(exp(-1i*alpha)*ZETA +(exp(1i*alpha)*r^2)./ZETA) + (1i*G*log(ZETA))/(2*pi);
    psiZeta = imag(FZeta);
    
    % Complex Velocity in the Zeta Domain 
    wZeta = U*(exp(-1i*alpha) - (exp(1i*alpha)*r^2)./(ZETA.*ZETA)) + (1i*G)./(2*pi*ZETA);
    uZeta = real(wZeta);
    vZeta = -imag(wZeta);
    Umag_Zeta = sqrt(wZeta.*conj(wZeta));
    
    % Plot the Velocity Field in Zeta
    fig = figure('Visible','off');
    subplot(2,2,1);
    hold on;
    grid on;
    contLevels = linspace(0,max(uZeta(:)),20);
    contourf(real(ZETA),imag(ZETA),Umag_Zeta,contLevels);
    colorbar;
    title("Velocity Contour in the $$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(2,2,2);
    hold on;
    grid on;
    contLevels = linspace(-1,1,20);
    contourf(real(ZETA),imag(ZETA),psiZeta,contLevels);
    colorbar;
    title("Streamfunction in the $$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(2,2,3);
    hold on;
    grid on;
    quiver(real(ZETA),imag(ZETA),uZeta,vZeta,'AutoScale','on','Color','b');
    title("Velocity Vectors in the $$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    
    hold off;
    outputFile = fullfile('./Results','09_ZetaDomainFlow.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Making a new figure
    magnitude = sqrt(uZeta.^2+vZeta.^2);
    uZetaNormalized = uZeta./magnitude;
    vZetaNormalized = vZeta./magnitude;
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(uZeta(:)),20);
    contourf(real(ZETA),imag(ZETA),Umag_Zeta,contLevels,'LineStyle','none');
    colorbar;
    hold on;
    quiver(real(ZETA),imag(ZETA),uZetaNormalized,vZetaNormalized,0.025,'k');
    title("Velocity Field - $$\zeta$$ Domain","Interpreter","latex");
    hold off;
    axis equal;

    outputFile = fullfile('./Results','10_VelocityField_ZetaDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    
    disp("Step 8A: Flow of the Zeta (Perfect Cirlce) Domain is computed.");
    
    %% STEP 8B - Velocity Field Computations in the Zeta Prime Domain
    % Refer Derivations.pdf 
    f1_zeta = 0;
    f2_zeta = 0;
    
    % Calculating auxillary functions 
    for m =1:length(A)
        f1_zeta = f1_zeta + (A(m)+1i*B(m))./(ZETA.^m);
        f2_zeta = f1_zeta + (m*(A(m)+1i*B(m)))./(ZETA.^m);
    end
    
    dZetaPrimedZeta = exp(f1_zeta).*(1-f2_zeta);
    wZetaPrime = wZeta./dZetaPrimedZeta;
    uZetaPrime = real(wZetaPrime);
    vZetaPrime = -imag(wZetaPrime);
    Umag_ZetaPrime = sqrt(wZetaPrime.*conj(wZetaPrime));
    
    % Plot the Velocity Field in Zeta Prime
    fig = figure('Visible','off');
    subplot(1,2,1);
    hold on;
    grid on;
    contLevels = linspace(0,max(uZetaPrime(:)),20);
    % The velocity field is computed on the centered zeta prime domain
    contourf(real(ZETA_PRIME),imag(ZETA_PRIME),Umag_ZetaPrime,contLevels);
    colorbar;
    title("Velocity Contour in the $$\zeta '$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(1,2,2);
    hold on;
    grid on;
    quiver(real(ZETA_PRIME),imag(ZETA_PRIME),uZetaPrime,vZetaPrime,'AutoScale','on','Color','b');
    title("Velocity Vectors in the $$\zeta '$$ Domain","Interpreter","latex");
    axis equal;
    
    hold off;

    outputFile = fullfile('./Results','11_ZetaPrimeDomainFlow.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Making a new figure
    magnitude = sqrt(uZetaPrime.^2+vZetaPrime.^2);
    uZetaPrimeNormalized = uZetaPrime./magnitude;
    vZetaPrimeNormalized = vZetaPrime./magnitude;
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(uZetaPrime(:)),20);
    contourf(real(ZETA_PRIME),imag(ZETA_PRIME),Umag_ZetaPrime,contLevels,'LineStyle','none');
    colorbar;
    hold on;
    quiver(real(ZETA_PRIME),imag(ZETA_PRIME),uZetaPrimeNormalized,vZetaPrimeNormalized,0.025,'k');
    title("Velocity Field - $$\zeta'$$ Domain","Interpreter","latex");
    hold off;
    axis equal;
    outputFile = fullfile('./Results','12_VelocityField_ZetaPrimeDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 8B: Flow of the Zeta Prime (Intermediate Cirlce) Domain is computed.");
    
    %% STEP 8C - Velocity Field Computations in the Z Domain
    % Compute the velocity field in the Z domain '
    dZdZetaPrime = 1 - a^2./(ZETA_PRIME.^2);
    wZ = wZeta./dZdZetaPrime;
    uZ = real(wZ);
    vZ = -imag(wZ);
    Umag_Z = sqrt(wZ.*conj(wZ));
    
    
    % Plot the Velocity Field in Zeta Prime
    fig = figure('Visible','off');
    subplot(1,2,1);
    hold on;
    grid on;
    contLevels = linspace(0,max(uZ(:)),20);
    contourf(real(Z),imag(Z),Umag_Z,contLevels);
    colorbar;
    title("Velocity Contour in the $$z$$ Domain","Interpreter","latex");
    axis equal;
    
    subplot(1,2,2);
    hold on;
    grid on;
    quiver(real(Z),imag(Z),uZ,vZ,'AutoScale','on','Color','b');
    title("Velocity Vectors in the $$z$$ Domain","Interpreter","latex");
    axis equal;
    
    hold off;
    outputFile = fullfile('./Results','13_ZDomainFlow.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);

    
    % Making a new figure
    magnitude = sqrt(uZ.^2+vZ.^2);
    uZNormalized = uZ./magnitude;
    vZNormalized = vZ./magnitude;
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(Umag_Z(:)),100);
    contourf(real(Z),imag(Z),Umag_Z,contLevels,'LineStyle','none');
    colorbar;
    hold on;
    quiver(real(Z),imag(Z),uZNormalized,vZNormalized,0.025,'k');
    title("Velocity Field - $$\mathcal{Z}$$ Domain","Interpreter","latex");
    hold off;
    axis equal;
    outputFile = fullfile('./Results','14_VelocityField_ZDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 8C: Flow of the Z (Perfect Cirlce) Domain is computed.");
    
    %% STEP 9 - Computing the Acceleration and Appellian in the Z and Zeta Domain
    %% STEP 9A - Computing the Acceleration and Appellian in Zeta Domain
    dwdZeta = (2*U*exp(1i*alpha)*r^2)./(ZETA.^3) + (1i*G)./(2*pi*(ZETA.^2));
    accZeta = wZeta.*conj(dwdZeta);
    accZetaMagSq = accZeta.*conj(accZeta);
    accZetaMag = sqrt(accZetaMagSq);
    
    % Plot the acceleration 
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(accZetaMag(:)),50);
    contourf(real(ZETA),imag(ZETA),accZetaMag,contLevels,'LineStyle','none');
    colorbar;
    title("Acceleration Magnitude in $$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','15_AccelerationMagnitude_ZetaDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Computing the grid and appellian by linear interpolation
    rows = radialPts - 1;
    columns = angularPts - 1;
    Xi = real(ZETA);
    Eta = imag(ZETA);
    Xic = zeros(rows,columns);
    Etac = zeros(rows,columns);
    accZetaMagSqCentroid = zeros(rows,columns);
    areaZetaC = zeros(rows,columns);
    for i = 1:rows
        for j =1:columns
            xi = [Xi(i,j) Xi(i,j+1) Xi(i+1,j+1) Xi(i+1,j)];
            eta = [Eta(i,j) Eta(i,j+1) Eta(i+1,j+1) Eta(i+1,j)];
            Xic(i,j) = sum(xi)/4;
            Etac(i,j) = sum(eta)/4;
            area = 0.5 * abs(xi(1)*eta(2) + xi(2)*eta(3) + xi(3)*eta(4) ...
                + xi(4)*eta(1) - (eta(1)*xi(2) + eta(2)*xi(3) + eta(3)*xi(4) ...
                + eta(4)*xi(1)));
            areaZetaC(i,j) = area;
            accZetaMagSqCentroid(i,j) = (accZetaMagSq(i,j)+accZetaMagSq(i,j+1) ...
                +accZetaMagSq(i+1,j+1)+accZetaMagSq(i+1,j))/4;
        end
    end
    
    % Assuming Density is 1
    appellian = 0.5*(accZetaMagSqCentroid.*areaZetaC);
    appellian = sum(appellian,"all");
    normApp = appellian/U^4;
    disp("The normalized appellian (for equivalent Kutta circulation) in the Zeta domain is given as: ");
    disp(normApp)
    
    %% STEP 9B - Computing the Acceleration and Appellian in the Z Domain
    
    % Computing d2ZdZetaPrime2
    d2ZdZetaPrime2 = (2*a^2)./(ZETA_PRIME.^3);
    
    % Compute d2zeta'dzeta2
    f3_zeta = 0;
    f4_zeta = 0; 
    for m =1:length(A)
        f3_zeta = f3_zeta + ((-m)*(A(m)+1i*B(m)))./(ZETA.^(m+1));
        f4_zeta = f4_zeta + (m*(m+1)*(A(m)+1i*B(m)))./(ZETA.^(m+2));
    end
    d2ZetaPrimedZeta2 = exp(f1_zeta).*(f3_zeta.*(1-f2_zeta)+f4_zeta);
    
    % Compute G_zeta, conjugate, dGdZeta
    G_zeta = 1./(dZdZetaPrime.*dZetaPrimedZeta);
    dZdZeta = dZdZetaPrime.*dZetaPrimedZeta;
    d2ZdZeta2 = d2ZdZetaPrime2.*(dZetaPrimedZeta.^2) + dZdZetaPrime.*d2ZetaPrimedZeta2;
    dGdZeta = -d2ZdZeta2./(dZdZeta.^2);
    
    % Compute the acceleration in Z domain
    accZMag = abs(G_zeta).*abs(conj(G_zeta).*accZeta + (abs(wZeta).^2).*conj(dGdZeta));
    accZMagSq = accZMag.^2;
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(accZMagSq(:)),100);
    contourf(real(Z),imag(Z),accZMagSq,contLevels);
    colorbar;
    title("Acceleration Squared Contour-Map in the $$z$$ Domain","Interpreter","latex");
    axis equal;
    outputFile = fullfile('./Results','16_AccelerationMagnitude_ZDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Computing the grid and appellian by linear interpolation
    rows = radialPts - 1;
    columns = angularPts - 1;
    X = real(Z);
    Y = imag(Z);
    Xc = zeros(rows,columns);
    Yc = zeros(rows,columns);
    accZMagSqCentroid = zeros(rows,columns);
    areaZC = zeros(rows,columns);
    for i = 1:rows
        for j =1:columns
            x = [X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)];
            y = [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)];
            Xc(i,j) = sum(x)/4;
            Yc(i,j) = sum(y)/4;
            area = 0.5 * abs(x(1)*y(2) + x(2)*y(3) + x(3)*y(4) ...
                + x(4)*y(1) - (y(1)*x(2) + y(2)*x(3) + y(3)*x(4) ...
                + y(4)*x(1)));
            areaZC(i,j) = area;
            accZMagSqCentroid(i,j) = (accZMagSq(i,j)+accZMagSq(i,j+1) ...
                +accZMagSq(i+1,j+1)+accZMagSq(i+1,j))/4;
        end
    end
    
    % Assuming Density is 1
    appellian = 0.5*(accZMagSqCentroid.*areaZC);
    appellian = sum(appellian,"all");
    normApp = appellian/U^4;
    disp("The normalized appellian (for equivalent Kutta circulation) in the Z domain is given as: ");
    disp(normApp)
    
    
    
    
    %% STEP 10: Variation of Appellian with Gamma
    
    % Initialize the Normalized Gamma
    kutta_normG = sin(alpha+beta); 
    normG = linspace(-kuttaNormGDev + kutta_normG,kuttaNormGDev + kutta_normG,101); %Always keep number of points as odd
    % Initialize a vector to store the Normalized appellian
    normAppVec = zeros(1,length(normG));
    
    
    for p = 1:length(normG)
        G = 4*pi*U*r*normG(p);
    
        % STEP 8A - Velocity Field Computations in the Zeta Domain
        % Complex Potential in the Zeta Domain
        FZeta = U*(exp(-1i*alpha)*ZETA +(exp(1i*alpha)*r^2)./ZETA) + (1i*G*log(ZETA))/(2*pi);
        psiZeta = imag(FZeta);    
        % Complex Velocity in the Zeta Domain 
        wZeta = U*(exp(-1i*alpha) - (exp(1i*alpha)*r^2)./(ZETA.*ZETA)) + (1i*G)./(2*pi*ZETA);
        uZeta = real(wZeta);
        vZeta = -imag(wZeta);
        Umag_Zeta = sqrt(wZeta.*conj(wZeta));
        
        % STEP 8B - Velocity Field Computations in the Zeta Prime Domain
        % Refer Derivations.pdf 
        f1_zeta = 0;
        f2_zeta = 0;
        for m =1:length(A)
            f1_zeta = f1_zeta + (A(m)+1i*B(m))./(ZETA.^m);
            f2_zeta = f1_zeta + (m*(A(m)+1i*B(m)))./(ZETA.^m);
        end
        dZetaPrimedZeta = exp(f1_zeta).*(1-f2_zeta);
        wZetaPrime = wZeta./dZetaPrimedZeta;
        uZetaPrime = real(wZetaPrime);
        vZetaPrime = -imag(wZetaPrime);
        Umag_ZetaPrime = sqrt(wZetaPrime.*conj(wZetaPrime));
        
        % STEP 8C - Velocity Field Computations in the Z Domain
        % Compute the velocity field in the Z domain '
        dZdZetaPrime = 1 - a^2./(ZETA_PRIME.^2);
        wZ = wZeta./dZdZetaPrime;
        uZ = real(wZ);
        vZ = -imag(wZ);
        Umag_Z = sqrt(wZ.*conj(wZ));
        
        % STEP 9B - Computing the Acceleration and Appellian in the Z Domain
        % Follow Derivations.pdf; Note: dZdZetaPrime and dZetaPrimedZeta are computed above
        
        % Computing d2ZdZetaPrime2
        d2ZdZetaPrime2 = (2*a^2)./(ZETA_PRIME.^3);
        
        % Compute d2zeta'dzeta2
        f3_zeta = 0;
        f4_zeta = 0; 
        for m =1:length(A)
            f3_zeta = f3_zeta + ((-m)*(A(m)+1i*B(m)))./(ZETA.^(m+1));
            f4_zeta = f4_zeta + (m*(m+1)*(A(m)+1i*B(m)))./(ZETA.^(m+2));
        end
        d2ZetaPrimedZeta2 = exp(f1_zeta).*(f3_zeta.*(1-f2_zeta)+f4_zeta);
        
        % Compute G_zeta, conjugate, dGdZeta
        G_zeta = 1./(dZdZetaPrime.*dZetaPrimedZeta);
        dZdZeta = dZdZetaPrime.*dZetaPrimedZeta;
        d2ZdZeta2 = d2ZdZetaPrime2.*(dZetaPrimedZeta.^2) + dZdZetaPrime.*d2ZetaPrimedZeta2;
        dGdZeta = -d2ZdZeta2./(dZdZeta.^2);
        
        % Compute the acceleration in Z domain
        accZMag = abs(G_zeta).*abs(conj(G_zeta).*accZeta + (abs(wZeta).^2).*conj(dGdZeta));
        accZMagSq = accZMag.^2;
        
        % Computing the grid and appellian by linear interpolation
        rows = radialPts - 1;
        columns = angularPts - 1;
        X = real(Z);
        Y = imag(Z);
        Xc = zeros(rows,columns);
        Yc = zeros(rows,columns);
        accZMagSqCentroid = zeros(rows,columns);
        areaZC = zeros(rows,columns);
        for i = 1:rows
            for j =1:columns
                x = [X(i,j) X(i,j+1) X(i+1,j+1) X(i+1,j)];
                y = [Y(i,j) Y(i,j+1) Y(i+1,j+1) Y(i+1,j)];
                Xc(i,j) = sum(x)/4;
                Yc(i,j) = sum(y)/4;
                area = 0.5 * abs(x(1)*y(2) + x(2)*y(3) + x(3)*y(4) ...
                    + x(4)*y(1) - (y(1)*x(2) + y(2)*x(3) + y(3)*x(4) ...
                    + y(4)*x(1)));
                areaZC(i,j) = area;
                accZMagSqCentroid(i,j) = (accZMagSq(i,j)+accZMagSq(i,j+1) ...
                    +accZMagSq(i+1,j+1)+accZMagSq(i+1,j))/4;
            end
        end
        
        % Assuming Density is 1
        appellian = 0.5*(accZMagSqCentroid.*areaZC);
        appellian = sum(appellian,"all");
        normApp = appellian/U^4;
    
        normAppVec(p) = normApp;
       
    end
    
    %% STEP 11: PLOT MINIMIZATION
    % Find the minimizing Gamma
    [minApp, idx] = min(normAppVec);
    minNormG = normG(idx);
    kuttaApp = normAppVec(uint8(length(normG)/2));
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    plot(normG,normAppVec,'b','LineWidth',2);
    xlabel("Normalized Circulation");
    ylabel("Normalized Appellian");
    title("Application of PMPG");
    hold on;
    scatter(minNormG, minApp, 'r', 'filled');
    scatter(kutta_normG, kuttaApp, 'k', 'filled');
    legend("Variation of Appellian","PMPG","Kutta Condition");
    hold off;
    outputFile = fullfile('./Results','17_Application_of_PMPG.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Kutta Circulation:");
    disp(kutta_normG);
    disp("Minimized Circulation by PMPG:");
    disp(minNormG);
    
    %% STEP 12 : PLOT THE FLOW FIELD WITH MINIMIZED CIRCULATION
    G = 4*pi*U*r*minNormG;
    wZeta = U*(exp(-1i*alpha) - (exp(1i*alpha)*r^2)./(ZETA.*ZETA)) + (1i*G)./(2*pi*ZETA);
    wZetaPrime = wZeta./dZetaPrimedZeta;
    wZ = wZetaPrime./dZdZetaPrime;
    uZ = real(wZ);
    vZ = -imag(wZ);
    Umag_Z = sqrt(wZ.*conj(wZ));
    magnitude = sqrt(uZ.^2+vZ.^2);
    uZNormalized = uZ./magnitude;
    vZNormalized = vZ./magnitude;
    
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(Umag_Z(:)),100);
    contourf(real(Z),imag(Z),Umag_Z,contLevels,'LineStyle','none');
    colorbar;
    hold on;
    quiver(real(Z),imag(Z),uZNormalized,vZNormalized,0.025,'k');
    title("Velocity Field - $$\mathcal{Z}$$ Domain","Interpreter","latex");
    hold off;
    axis equal;
    outputFile = fullfile('./Results','18_VelocityField_PMPG_Solution.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);

    flowFieldData.minNormG = minNormG;
    flowFieldData.minApp = minApp;
    flowFieldData.kuttaApp = kuttaApp;
    flowFieldData.kuttaNormG = kutta_normG;
end