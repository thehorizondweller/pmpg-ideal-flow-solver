function [x1, x2, x3, h, tl] = perSpline(theta, rho)
% Revised for star-like parametrization rho = rho(theta) for Theodorsen
    x = rho(:);
    if abs(x(1) - x(end)) > 100*eps
        x(end+1) = x(1);
    end
    n1 = length(x);
    n = n1 - 1;
    dx = diff(x);
    h = diff(theta);
    h = h(:);
    tl = sum(h);
    h(n1) = h(1);
    p = h(1:n);
    q = h(2:n1);
    a = q./(p + q);
    b = 1 - a;
    c = spdiags([[b(n); ones(n-1,1)] [a(2:n);0] [2*ones(n,1)] [0;b(1:n-1)] ...
        [ones(n-1,1); a(1)]], ...
        [-n+1 -1 0 1 n-1], n, n);
    d1 = 3*(a.*dx./p + b.*[dx(2:n); x(2) - x(n1)]./q);
    mmdflag = spparms('autommd');
    spparms('autommd', 0);
    x1 = c\d1;
    spparms('autommd', mmdflag);
    x1(2:n1) = x1;
    x1(1) = x1(n1);
    x2(2:n1) = 2*(x1(1:n) + 2*x1(2:n1) - 3*dx./p)./p;
    x2(1) = x2(n1);
    x2 = x2';
    x3 = diff(x2)./p;
    x3(n1) = x3(1);
end
