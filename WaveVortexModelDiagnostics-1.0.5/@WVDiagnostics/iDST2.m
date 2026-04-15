function matrix = iDST2(N)
% InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix.
%
% InverseSineTransformMatrix_DST2  Inverse of forward DST-II matrix
% This inverts SineTransformForwardMatrix_DST2 when the forward uses the
% 2/N scaling above (with no endpoint tweaks).
% Implementation: DST-III with the LAST column halved (k = N).
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = iDST2(N)
% - Parameter N: input argument `N`
% - Returns matrix: output value `matrix`
arguments
    N
end


matrix = zeros(N,N);

for j = 1:N
    for k = 1:N
        matrix(j,k) = sin(pi*(j-0.5)*k/N);
    end
end

% Halve the highest-frequency column (k = N)
matrix(:,N) = matrix(:,N)/2;

return