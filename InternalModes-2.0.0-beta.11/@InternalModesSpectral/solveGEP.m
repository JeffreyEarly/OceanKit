function [V_cheb,h] = solveGEP(self,A,B,options)
% Solve and sort a generalized eigenvalue problem.
%
% - Topic: Developer topics
% - Developer: true
% - Declaration: [V_cheb,h] = solveGEP(self,A,B,options)
% - Parameter self: InternalModesSpectral instance
% - Parameter A: left generalized-eigenproblem matrix
% - Parameter B: right generalized-eigenproblem matrix
% - Parameter options.negativeEigenvalues: number of negative-depth modes to retain first
% - Returns V_cheb: sorted eigenvectors in Chebyshev coefficient space
% - Returns h: sorted equivalent-depth vector
arguments
    self InternalModesSpectral
    A (:,:) double
    B (:,:) double
    options.negativeEigenvalues (1,1) double = 0
end

if any(isnan(A(:))) || any(isnan(B(:)))
    error('InternalModesSpectral:NaNInMatrix', 'EVP setup failed. Found at least one NaN in matrices A and B.');
end

[V,D] = eig(A, B);
[h, permutation] = sort(real(self.hFromLambda(diag(D))), 'descend');
V_cheb = V(:,permutation);

if options.negativeEigenvalues > 0
    negIndices = find(h < 0, options.negativeEigenvalues, 'first');
    permutation = cat(1, negIndices, setdiff((1:length(h))', negIndices));
    h = h(permutation);
    V_cheb = V_cheb(:,permutation);
end
end
