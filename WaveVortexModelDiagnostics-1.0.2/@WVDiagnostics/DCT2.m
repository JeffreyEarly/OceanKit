function matrix = DCT2(N)
% CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix.
%
% CosineTransformForwardMatrix_DCT2  Discrete Cosine Transform (DCT-II) matrix
% Forward scaling: 2/N (no extra endpoint halving here).
% With this choice, the inverse (DCT-III) is the plain cosine matrix with
% the DC (first) COLUMN halved.
% If X = CosineTransformForwardMatrix_DCT2(N) * x, then
% x = InverseCosineTransformMatrix_DCT2(N) * X.
% See also: InverseCosineTransformMatrix_DCT2
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = DCT2(N)
% - Parameter N: input argument `N`
% - Returns matrix: output value `matrix`
arguments
    N
end


matrix = zeros(N,N);

for k = 1:N
    for j = 1:N
        matrix(k,j) = (2/N) * cos(pi*(j-0.5)*(k-1)/N);
    end
end

return