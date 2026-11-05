function [z,weights] = neutralGeostrophicQuadrature(joins,nPoints)
% Return independent Clenshaw--Curtis quadrature on each physical interval.
arguments
    joins (1,:) double {mustBeReal,mustBeFinite}
    nPoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nPoints,3)}
end
if any(diff(joins) <= 0)
    error("NeutralGeostrophicReference:InvalidJoins","Physical interval joins must increase strictly.");
end
n = nPoints-1;
theta = (0:n)'*pi/n;
x = cos(theta);
w = zeros(n+1,1);
interior = (2:n)';
v = ones(n-1,1);
if mod(n,2) == 0
    w([1 end]) = 1/(n*n-1);
    for k = 1:n/2-1
        v = v-2*cos(2*k*theta(interior))/(4*k*k-1);
    end
    v = v-cos(n*theta(interior))/(n*n-1);
else
    w([1 end]) = 1/(n*n);
    for k = 1:(n-1)/2
        v = v-2*cos(2*k*theta(interior))/(4*k*k-1);
    end
end
w(interior) = 2*v/n;
z = zeros(0,1);
weights = zeros(0,1);
for i = 1:numel(joins)-1
    interval = joins(i:i+1);
    zi = flipud(mean(interval)+diff(interval)*x/2);
    wi = flipud(diff(interval)*w/2);
    if i == 1
        z = zi;
        weights = wi;
    else
        weights(end) = weights(end)+wi(1);
        z = [z;zi(2:end)]; %#ok<AGROW>
        weights = [weights;wi(2:end)]; %#ok<AGROW>
    end
end
end
