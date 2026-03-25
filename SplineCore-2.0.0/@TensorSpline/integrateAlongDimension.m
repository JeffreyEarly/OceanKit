function [xi, tKnot, K] = integrateAlongDimension(xi, tKnot, K, dim, xMean, xStd)
% Integrate tensor coefficients along one dimension.
%
% - Topic: Utility
% - Developer: true
% - Declaration: [xi,tKnot,K] = integrateAlongDimension(xi,tKnot,K,dim,xMean,xStd)
perm = [dim, 1:(dim-1), (dim+1):ndims(xi)];
xiPermuted = permute(xi, perm);
xiMatrix = reshape(xiPermuted, size(xiPermuted, 1), []);

originalK = K;
originalTKnot = tKnot;
transformedSlice = cumsum(BSpline(S=originalK-1, knotPoints=originalTKnot, xi=xiMatrix(:,1), xMean=xMean, xStd=xStd));
tKnot = transformedSlice.knotPoints;
K = transformedSlice.K;

transformedMatrix = zeros(numel(transformedSlice.xi), size(xiMatrix, 2), 'like', transformedSlice.xi);
transformedMatrix(:,1) = transformedSlice.xi;
for iSlice = 2:size(xiMatrix, 2)
    transformedSlice = cumsum(BSpline(S=originalK-1, knotPoints=originalTKnot, xi=xiMatrix(:,iSlice), xMean=xMean, xStd=xStd));
    transformedMatrix(:,iSlice) = transformedSlice.xi;
end

outputSize = size(xiPermuted);
outputSize(1) = size(transformedMatrix, 1);
xi = ipermute(reshape(transformedMatrix, outputSize), perm);
