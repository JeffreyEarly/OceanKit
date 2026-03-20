function matrix = DST1(N)
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
%
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix
% This matrix exactly matches CosineTransformForward. See its documentation
% for details.
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = DST1(N)
% - Parameter N: input argument `N`
% - Returns matrix: output value `matrix`
arguments
    N
end


matrix = zeros(N-2,N);

for k=2:(N-1)
	for j=1:N
		matrix(k-1,j) = (2/(N-1))*sin(pi*(j-1)*(k-1)/(N-1));	
	end
end

return