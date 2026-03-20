function matrix = DST2(N)
% SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix.
%
% SineTransformForwardMatrix_DST2  Discrete Sine Transform (DST-II) matrix
% Forward scaling: 2/N (no endpoint tweaks).
% With this choice, the inverse (DST-III) is the plain sine matrix with
% the HIGHEST-frequency (last) COLUMN halved.
% If X = SineTransformForwardMatrix_DST2(N) * x, then
% x = InverseSineTransformMatrix_DST2(N) * X.
% See also: InverseSineTransformMatrix_DST2
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = DST2(N)
% - Parameter N: input argument `N`
% - Returns matrix: output value `matrix`
arguments
    N
end


matrix = zeros(N,N);

for k = 1:N
    for j = 1:N
        matrix(k,j) = (2/N) * sin(pi*(j-0.5)*k/N);
    end
end

return