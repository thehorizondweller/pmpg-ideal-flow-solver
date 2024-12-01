function flowFieldData = computeKTFlowField(radialPts,angularPts,terminalR, ...
    a,psi_0,zetaPrimeOffset,zetaPrime1,zetaPrime2,z1,z2,U,alpha,beta,kuttaNormGDev,A,B,n,zRetrace,zeta)
    
    %% STEP 7 - Build the Domains
    
    % Build the Domain

    % USER INPUTS - Manually update only for debugging
    % radialPts = 101;
    % angularPts = 501;
    % terminalR = 10;

    r = a*exp(psi_0); % Radius of the circle in the zeta domain
    domainGrids = flowDomainBuild_KTMap(radialPts, angularPts, terminalR, ...
        r, A, B, zetaPrimeOffset, zetaPrime1, zetaPrime2, z1, z2, beta); %, offset);
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

    outputFile = fullfile('./Results','09_Grid_Domains.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    disp("Step 7: Grid Constructed for all 3 domains.")
    
    %% STEP 8A - Velocity Field Computations in the Zeta Domain

    % USER INPUTS - Manually update for debugging only 
    % U = 1; % Far Field Conditions
    % alpha = 5; % In degrees

    alpha = alpha * pi / 180; % Convert to Radians
    
    % Find the TE point 
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
        l1 = zRetrace(prev) - zRetrace(i);
        l2 = zRetrace(next) - zRetrace(i);
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
    outputFile = fullfile('./Results','10_Variation_of_Turning_Angle_ComputedDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    % Initial Guess for Gamma based on the Kutta Condition
    [minAngle, minIdx] = min(angles);
    TEAngleZetaDomain = angle(zeta(minIdx));
    kutta_normG = sin(alpha-TEAngleZetaDomain);
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

    outputFile = fullfile('./Results','11_KuttaCondFlowField_ZetaDomain.png');
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
    outputFile = fullfile('./Results','12_KuttaCond_VelocityField_ZetaDomain.png');
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

    outputFile = fullfile('./Results','13_KuttaCond_FlowField_ZetaPrimeDomain.png');
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

    outputFile = fullfile('./Results','14_KuttaCond_VelocityField_ZetaPrimeDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
    
    disp("Step 8B: Flow of the Zeta Prime (Intermediate Cirlce) Domain is computed.");
    
    %% STEP 8C - Velocity Field Computations in the Z Domain
    % Compute the velocity field in the Z domain
    q = ((ZETA_PRIME-zetaPrime1)./(ZETA_PRIME-zetaPrime2)).^beta;
    
    % CAUTION: Branch monitoring is required to compute
    % (zeta'-zeta'_1)^(beta-1) which we call root term 1
    [rowNum, colNum] = size(ZETA_PRIME);
    rootTerm1 = zeros(rowNum, colNum);
    for i=1:rowNum
        rootTerm1(i,1) = (ZETA_PRIME(i,1)-zetaPrime1).^(beta-1);
        rootTerm1(i,2) = (ZETA_PRIME(i,2)-zetaPrime1).^(beta-1);
        w = rootTerm1(i,2);
        for j=3:colNum
            branch = w;
            w = ZETA_PRIME(i,j)-zetaPrime1;
            w = crootbr(w, branch, beta-1);
            rootTerm1(i,j) = w;
        end
    end
    
    dqdZetaPrime = beta*(zetaPrime1-zetaPrime2)*rootTerm1./((ZETA_PRIME-zetaPrime2).^(beta+1));
    dZdZetaPrime = (z1-z2)*dqdZetaPrime./((1-q).^2);
    
    % Scale dZdZetaPrime to go to 1 as z or zeta' goes to inf
    sf = (1/beta)*((z1-z2)./(zetaPrime1-zetaPrime2));
    dZdZetaPrime_scaled = dZdZetaPrime*(1/sf);
    wZ = wZetaPrime./dZdZetaPrime_scaled;
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
    hold on;
    markerSize = 10; % Size of markers
    markerColor = 'r'; % Color of markers (e.g., 'r' for red)
    markerShape = 's'; % Shape of markers (e.g., 's' for square)
    scatter(real(zRetrace),imag(zRetrace),markerSize,markerColor,'filled',markerShape);
    title("Velocity Vectors in the $$z$$ Domain","Interpreter","latex");
    axis equal;
    
    hold off;

    outputFile = fullfile('./Results','15_KuttaCond_FlowField_ZDomain.png');
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
    
    outputFile = fullfile('./Results','16_KuttaCond_VelocityField_ZDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);

    disp("Step 8C: Flow of the Z (Perfect Cirlce) Domain is computed.");
    
    %% STEP 9A: Computing Acceleration and Appellian in the ZETA Domain
    dwdZeta = (2*U*exp(1i*alpha)*r^2)./(ZETA.^3) + (1i*G)./(2*pi*(ZETA.^2));
    accZeta = wZeta.*conj(dwdZeta);
    accZetaMagSq = accZeta.*conj(accZeta);
    accZetaMag = sqrt(accZetaMagSq);
    
    % Plot the acceleration 
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(accZetaMag(:)),20);
    contourf(real(ZETA),imag(ZETA),accZetaMag,contLevels,'LineStyle','none');
    colorbar;
    title("Acceleration Magnitude in $$\zeta$$ Domain","Interpreter","latex");
    axis equal;
    hold off;

    outputFile = fullfile('./Results','17_KuttaCond_AccMag_ZetaDomain.png');
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
    disp("The normalized appellian (for Kutta condition) in the Zeta domain is given as: ");
    disp(normApp)
    
    %% STEP 9b: Computing the Acceleration and Appellian in the Z domain
    % Compute d2zeta'dzeta2
    f3_zeta = 0;
    f4_zeta = 0; 
    for m =1:length(A)
        f3_zeta = f3_zeta + ((-m)*(A(m)+1i*B(m)))./(ZETA.^(m+1));
        f4_zeta = f4_zeta + (m*(m+1)*(A(m)+1i*B(m)))./(ZETA.^(m+2));
    end
    d2ZetaPrimedZeta2 = exp(f1_zeta).*(f3_zeta.*(1-f2_zeta)+f4_zeta);
    
    % Compute d2ZdZetaPrime2
    % CAUTION: Branch monitoring is required to compute
    % (zeta'-zeta'_1)^(beta-2) which we call root term 2
    [rowNum, colNum] = size(ZETA_PRIME);
    rootTerm2 = zeros(rowNum, colNum);
    for i=1:rowNum
        rootTerm2(i,1) = (ZETA_PRIME(i,1)-zetaPrime1).^(beta-2);
        rootTerm2(i,2) = (ZETA_PRIME(i,2)-zetaPrime1).^(beta-2);
        w = rootTerm2(i,2);
        for j=3:colNum
            branch = w;
            w = ZETA_PRIME(i,j)-zetaPrime1;
            w = crootbr(w, branch, beta-2);
            rootTerm2(i,j) = w;
        end
    end
    d2qdZetaPrime2 = (beta*(zetaPrime1-zetaPrime2)*rootTerm2./((ZETA_PRIME-zetaPrime2).^(beta+2))).*(-ZETA_PRIME + beta*zetaPrime1 - (beta-1)*zetaPrime2);
    d2ZdZetaPrime2 = ((z1-z2)/sf)*((-2*dqdZetaPrime.^2)./((1-q).^3) + d2qdZetaPrime2./((1-q).^2));
    
    % Compute G(zeta), conjugate and dGdZeta
    G_zeta = 1./(dZdZetaPrime_scaled.*dZetaPrimedZeta);
    dZdZeta = dZdZetaPrime_scaled.*dZetaPrimedZeta;
    d2ZdZeta2 = d2ZdZetaPrime2.*(dZetaPrimedZeta.^2) + dZdZetaPrime_scaled.*d2ZetaPrimedZeta2;
    dGdZeta = -d2ZdZeta2./(dZdZeta.^2);
    
    % Compute the acceleration in Z domain
    accZMag = abs(G_zeta).*abs(conj(G_zeta).*accZeta + (abs(wZeta).^2).*conj(dGdZeta));
    accZMagSq = accZMag.^2;
    
    % Plot the acceleration 
    fig = figure('Visible','off');
    hold on;
    grid on;
    contLevels = linspace(0,max(accZMag(:)),20);
    contourf(real(Z),imag(Z),accZMag,contLevels,'LineStyle','none');
    colorbar;
    title("Acceleration Magnitude in $$\mathcal{Z}$$ Domain","Interpreter","latex");
    axis equal;
    hold off;
    outputFile = fullfile('./Results','18_KuttaCond_AccMag_ZDomain.png');
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
    disp("The normalized appellian (for the Kutta condition) in the Z domain is given as: ");
    disp(normApp)
    
    %% STEP 10: Minimization of Appellian 
    % Declare the FarField conditions
    % U = 1;
    % alpha = 5;
    % alpha = alpha * pi/180;
    
    % Initialize the Normalized Gamma
    % kutta_normG = sin(alpha-TEAngleZetaDomain);
    normG = linspace(-kuttaNormGDev+kutta_normG,kutta_normG+kuttaNormGDev,101); %Always keep number of points as odd
    
    % Initialize a vector to store the Normalized appellian
    normAppVec = zeros(1,length(normG));
    
    
    for p = 1:length(normG)
        G = 4*pi*U*r*normG(p);
    
        % STEP 8A - Velocity Field Computations in the Zeta Domain
        % Complex Velocity in the Zeta Domain 
        wZeta = U*(exp(-1i*alpha) - (exp(1i*alpha)*r^2)./(ZETA.*ZETA)) + (1i*G)./(2*pi*ZETA);
        uZeta = real(wZeta);
        vZeta = -imag(wZeta);
        Umag_Zeta = sqrt(wZeta.*conj(wZeta));
        
        % STEP 8B - Velocity Field Computations in the Zeta Prime Domain
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
        
        % STEP 8C - Velocity Field Computations in the Z Domain
        % Compute the velocity field in the Z domain
        q = ((ZETA_PRIME-zetaPrime1)./(ZETA_PRIME-zetaPrime2)).^beta;
        % CAUTION: Branch monitoring is required to compute
        % (zeta'-zeta'_1)^(beta-1) which we call root term 1
        [rowNum, colNum] = size(ZETA_PRIME);
        rootTerm1 = zeros(rowNum, colNum);
        for i=1:rowNum
            rootTerm1(i,1) = (ZETA_PRIME(i,1)-zetaPrime1).^(beta-1);
            rootTerm1(i,2) = (ZETA_PRIME(i,2)-zetaPrime1).^(beta-1);
            w = rootTerm1(i,2);
            for j=3:colNum
                branch = w;
                w = ZETA_PRIME(i,j)-zetaPrime1;
                w = crootbr(w, branch, beta-1);
                rootTerm1(i,j) = w;
            end
        end
        dqdZetaPrime = beta*(zetaPrime1-zetaPrime2)*rootTerm1./((ZETA_PRIME-zetaPrime2).^(beta+1));
        dZdZetaPrime = (z1-z2)*dqdZetaPrime./((1-q).^2);
        % Scale dZdZetaPrime to go to 1 as z or zeta' goes to inf
        sf = (1/beta)*((z1-z2)./(zetaPrime1-zetaPrime2));
        dZdZetaPrime_scaled = dZdZetaPrime*(1/sf);
        wZ = wZetaPrime./dZdZetaPrime_scaled;
        uZ = real(wZ);
        vZ = -imag(wZ);
        Umag_Z = sqrt(wZ.*conj(wZ));
        
        % STEP 9B: Computing the Acceleration and Appellian in the Z domain
        % Compute d2zeta'dzeta2
        f3_zeta = 0;
        f4_zeta = 0; 
        for m =1:length(A)
            f3_zeta = f3_zeta + ((-m)*(A(m)+1i*B(m)))./(ZETA.^(m+1));
            f4_zeta = f4_zeta + (m*(m+1)*(A(m)+1i*B(m)))./(ZETA.^(m+2));
        end
        d2ZetaPrimedZeta2 = exp(f1_zeta).*(f3_zeta.*(1-f2_zeta)+f4_zeta);
        
        % Compute d2ZdZetaPrime2
        % CAUTION: Branch monitoring is required to compute
        % (zeta'-zeta'_1)^(beta-2) which we call root term 2
        [rowNum, colNum] = size(ZETA_PRIME);
        rootTerm2 = zeros(rowNum, colNum);
        for i=1:rowNum
            rootTerm2(i,1) = (ZETA_PRIME(i,1)-zetaPrime1).^(beta-2);
            rootTerm2(i,2) = (ZETA_PRIME(i,2)-zetaPrime1).^(beta-2);
            w = rootTerm2(i,2);
            for j=3:colNum
                branch = w;
                w = ZETA_PRIME(i,j)-zetaPrime1;
                w = crootbr(w, branch, beta-2);
                rootTerm2(i,j) = w;
            end
        end
        d2qdZetaPrime2 = (beta*(zetaPrime1-zetaPrime2)*rootTerm2./((ZETA_PRIME-zetaPrime2).^(beta+2))).*(-ZETA_PRIME + beta*zetaPrime1 - (beta-1)*zetaPrime2);
        d2ZdZetaPrime2 = ((z1-z2)/sf)*((-2*dqdZetaPrime.^2)./((1-q).^3) + d2qdZetaPrime2./((1-q).^2));
        
        % Compute G(zeta), conjugate and dGdZeta
        G_zeta = 1./(dZdZetaPrime_scaled.*dZetaPrimedZeta);
        dZdZeta = dZdZetaPrime_scaled.*dZetaPrimedZeta;
        d2ZdZeta2 = d2ZdZetaPrime2.*(dZetaPrimedZeta.^2) + dZdZetaPrime_scaled.*d2ZetaPrimedZeta2;
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
    outputFile = fullfile('./Results','19_Application_of_PMPG.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);
    
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

    outputFile = fullfile('./Results','20_PMPG_VelocityField_ZDomain.png');
    exportgraphics(fig, outputFile, 'Resolution', 600);

    % Print the results
    fprintf("The normalized circulation which satisfies Kutta condition is: %4.2f \n", kutta_normG);
    fprintf("The normalized appellian (for Minimizing Circulation - PMPG) in the Z domain is given as: %4.2f \n", minApp);
    fprintf("The normalized circulation which minimizes the Appellian in the Z domain is: %4.2f \n", minNormG);

    flowFieldData.minNormG = minNormG;
    flowFieldData.minApp = minApp;
    flowFieldData.kuttaApp = kuttaApp;
    flowFieldData.kuttaNormG = kutta_normG;

end