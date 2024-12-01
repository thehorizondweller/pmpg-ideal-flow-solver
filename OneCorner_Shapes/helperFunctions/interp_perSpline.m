function fval = interp_perSpline(x, x1, x2, x3, h, tl, s)
% Revised for a single periodic spline for starlike parametrization
    x = x(:);
    s = s(:);
    n1 = length(x);
    nn = n1 - 1;
    n = length(s);
    s = s - floor(s/tl)*tl;
    cumsumh = [0; cumsum(h(1:nn))];
    ppx = mkpp(cumsumh, [x3(1:nn)/6; x2(1:nn)/2; x1(1:nn); x(1:nn)]);
    fval = ppval(ppx, s);
end
