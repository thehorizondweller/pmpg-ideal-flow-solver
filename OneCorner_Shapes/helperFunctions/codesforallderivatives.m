% ZetaPrime-Zeta
f1_zeta = 0;
f2_zeta = 0;
f3_zeta = 0;
f4_zeta = 0;

% Calculating auxillary functions 
for m =1:length(A)
    f1_zeta = f1_zeta + (A(m)+1i*B(m))./(ZETA.^m);
    f2_zeta = f1_zeta + (m*(A(m)+1i*B(m)))./(ZETA.^m);
    f3_zeta = f3_zeta + ((-m)*(A(m)+1i*B(m)))./(ZETA.^(m+1));
    f4_zeta = f4_zeta + (m*(m+1)*(A(m)+1i*B(m)))./(ZETA.^(m+2));
end

dZetaPrimedZeta = exp(f1_zeta).*(1-f2_zeta);
d2ZetaPrimedZeta2 = exp(f1_zeta).*(f3_zeta.*(1-f2_zeta)+f4_zeta);


% Z-ZetaPrime
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