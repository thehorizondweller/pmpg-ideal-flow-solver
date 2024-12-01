function airfoilGrid = kjSmoothTEAirfoilGrid(n, C, e, beta, D, radialPts, angularPts, terminalR)
%{
1. n = number of points on the circle
2. r = radius of the circle
3. C = scaling factor 
4. e = eccentricity
5. beta = trailing edge angle 
6. D = Smoothing Parameter
7. radialPts = Discretization in the Radial Direction
8. angularPts = Discretization in the Angular Direction
9. terminalR = Domain Size
%}
    r = C*(1+e);
    x_offset = e*C;
    y_offset = C*tan(beta);
    mu = -x_offset + 1i*y_offset;

    zeta = zeros(1,n);

    for i = 0:n-1
        theta = (i/n)*2*pi;
        zeta(i+1) = r*cos(theta) + 1i*r*sin(theta) + mu;
    end

    % Build the Grid
    radialLocs = logspace(log(r)/log(10), log(terminalR)/log(10), radialPts);
    angularLocs = linspace(0, 2*pi, angularPts);
    [R,T] = ndgrid(radialLocs,angularLocs); % generate the polar grid
    XI = R.*cos(T) + real(mu);
    ETA = R.*sin(T) + imag(mu);
    ZETA = XI + 1i*ETA; % Cartesian Grid of the Canonical Domain 

    smoothParam = (1-D)/(1+D);

    z = zeta + (smoothParam * C^2)./zeta; % Surface Points
    Z = ZETA + (smoothParam * C^2)./ZETA; % Grid Points

    % Output Data
    airfoilGrid.z = z;
    airfoilGrid.zeta = zeta;
    airfoilGrid.zetaCenter = mu;
    airfoilGrid.ZETA = ZETA;
    airfoilGrid.Z = Z;
end