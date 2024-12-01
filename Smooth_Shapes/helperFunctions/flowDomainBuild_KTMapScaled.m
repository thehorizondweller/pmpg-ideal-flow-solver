%{
flowDomainBuild.m - Function File
- Builds the zeta, zeta prime and z domain grids
- Takes as input the radial and angular discretization along with maximum
    domain radius
- Other inputs include parameters required for computing the conformal map
    i.e. the fourier coefficients and the offset 
%}

function domainGrids = flowDomainBuild_KTMapScaled(radialPts, angularPts, terminalR, r, A, B, zetaPrimeOffset, zetaPrime1, zetaPrime2, z1, z2, beta)%, offset)
    % Build the Zeta Domain
    n = radialPts;
    zeta = zeros(1,n);
    for i = 0:n-1
        theta = (i/n)*2*pi;
        zeta(i+1) = r*cos(theta) + 1i*r*sin(theta);
    end

    % Build the Grid
    radialLocs = logspace(log(r)/log(10), log(terminalR)/log(10), radialPts);
    angularLocs = linspace(0, 2*pi, angularPts);
    [R,T] = ndgrid(radialLocs,angularLocs); % generate the polar grid
    XI = R.*cos(T);
    ETA = R.*sin(T);
    ZETA = XI + 1i*ETA; % Cartesian Grid of the Canonical Domain
    
    % Build the Zeta Prime Domain
    zetaprime = zeros(1,n);
    for i = 1:n
        expSum = 0;
        for m = 1:length(A)
            expSum = expSum + (A(m)+1i*B(m))/zeta(i)^m;
        end
        zetaprime(i) = zeta(i)*exp(expSum);
    end

    ZETA_PRIME = ZETA;
    for i = 1:length(A)
        ZETA_PRIME = ZETA_PRIME.*exp((A(i)+1i*B(i))./(ZETA.^i));
    end
    
    % Adjust for the Offset - Uncentering
    ZETA_PRIME = ZETA_PRIME + zetaPrimeOffset;

    % Build the Z domain
    Z = inverseKTMapScaled(ZETA_PRIME, zetaPrime1, zetaPrime2, z1, z2, beta);

    domainGrids.ZETA = ZETA;
    domainGrids.ZETA_PRIME = ZETA_PRIME;
    domainGrids.Z = Z;

end