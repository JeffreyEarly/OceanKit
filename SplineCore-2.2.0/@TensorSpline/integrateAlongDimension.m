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
[transformedMatrix, tKnot, SIntegrated] = BSpline.integratedSplineState( ...
    xiMatrix, knotPoints=originalTKnot, S=originalK - 1, xMean=xMean, xStd=xStd);
K = SIntegrated + 1;

outputSize = size(xiPermuted);
outputSize(1) = size(transformedMatrix, 1);
xi = ipermute(reshape(transformedMatrix, outputSize), perm);
