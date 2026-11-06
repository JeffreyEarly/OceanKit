function [negativeCount,zeroPresent,matrix] = IMGeostrophicAPVEndpointInertia(IN,betaSurface,betaBottom,surfaceActive,bottomActive)
% Assess APV endpoint inertia without dividing by the stratification.
% IN is the integral of N2; beta values are reciprocal accelerations.
if surfaceActive && bottomActive
    matrix = [1+IN*betaSurface -1; -1 1+IN*betaBottom];
    % Compute the small eigenvalue from a factored determinant. Forming
    % 1+IN*beta first loses the finite external branch as IN tends to zero.
    sumBeta = betaSurface+betaBottom;
    productBeta = IN*betaSurface*betaBottom;
    determinantFactor = sumBeta+productBeta;
    zeroPresent = abs(determinantFactor)<=1e3*eps*(abs(sumBeta)+abs(productBeta));
    trace = 2+IN*sumBeta;
    radius = hypot(IN*(betaSurface-betaBottom)/2,1);
    direction = sign(trace);
    if direction==0, direction = 1; end
    large = trace/2+direction*radius;
    small = IN*determinantFactor/large;
    negativeCount = double(large<0)+double(~zeroPresent && small<0);
    return
elseif bottomActive
    matrix = 1+IN*betaBottom;
elseif surfaceActive
    matrix = 1+IN*betaSurface;
else
    matrix = zeros(0);
end
values = eig(matrix);
tolerance = 1e3*eps(max([1;abs(values(:))]));
negativeCount = nnz(values < -tolerance);
zeroPresent = any(abs(values) <= tolerance);
end
