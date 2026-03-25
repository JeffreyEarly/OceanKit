function spline = zeroSplineForDomain(domain, numDimensions, options)
% Create a zero spline over the supplied domain.
%
% - Topic: Utility
% - Developer: true
% - Declaration: spline = zeroSplineForDomain(domain,numDimensions,options)
arguments
    domain (:,2) double {mustBeNumeric,mustBeReal,mustBeFinite}
    numDimensions (1,1) double {mustBeInteger,mustBePositive}
    options.xStd = 1
end

tKnot = cell(1, numDimensions);
for iDim = 1:numDimensions
    tKnot{iDim} = reshape(domain(iDim,:), [], 1);
end
spline = TensorSpline(S=zeros(1, numDimensions), knotPoints=tKnot, xi=0, xMean=0, xStd=options.xStd);
