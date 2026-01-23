function matrix = iDCT1(n)
% CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix.
%
% CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix
% This matrix exactly matches CosineTransformBack. See its documentation
% for details.
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = iDCT1(n)
% - Parameter n: input argument `n`
% - Returns matrix: output value `matrix`
arguments
    n
end


matrix = zeros(n,n);

for j=1:n
	for k=1:n
		matrix(k,j) = cos(pi*(j-1)*(k-1)/(n-1));	
	end
end

matrix(:,1) = matrix(:,1)/2;

return