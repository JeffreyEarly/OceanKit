function matrix = DCT1(N)
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix.
%
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix
% This matrix exactly matches CosineTransformForward. See its documentation
% for details.
%
% - Topic: Transformations — Cosine and Sine
% - Declaration: matrix = DCT1(N)
% - Parameter N: input argument `N`
% - Returns matrix: output value `matrix`
arguments
    N
end


matrix = zeros(N,N);

for k=1:N
	for j=1:N
		matrix(k,j) = (2/(N-1))*cos(pi*(j-1)*(k-1)/(N-1));	
	end
end

matrix(1,:) = matrix(1,:)/2;
matrix(N,:) = matrix(N,:)/2;
matrix(:,1) = matrix(:,1)/2;
matrix(:,N) = matrix(:,N)/2;
matrix(1,:) = matrix(1,:)*2;

return