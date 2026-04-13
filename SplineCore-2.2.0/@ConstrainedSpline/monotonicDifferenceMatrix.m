function A = monotonicDifferenceMatrix(basisSize, dim, direction)
% Build coefficient-difference inequalities along one dimension.
%
% This helper constructs a sparse matrix `A` such that `A*xi <= 0`
% enforces adjacent coefficient differences of one sign along the selected
% tensor dimension.
%
% For the increasing case, each row encodes
%
% $$
% \xi_{j_1,\ldots,j_d} - \xi_{j_1,\ldots,j_k+1,\ldots,j_d} \le 0,
% $$
%
% so adjacent coefficients are nondecreasing along dimension `dim`. The
% decreasing case flips the sign pattern.
%
% - Topic: Compile constraints
% - Developer: true
% - Declaration: A = monotonicDifferenceMatrix(basisSize,dim,direction)
% - Parameter basisSize: tensor basis size per dimension
% - Parameter dim: constrained tensor dimension
% - Parameter direction: "increasing" or "decreasing"
% - Returns A: sparse inequality matrix acting on xi(:)
basisSize = reshape(basisSize, 1, []);
numCoefficients = prod(basisSize);
if isscalar(basisSize)
    basisSize = [basisSize, 1];
end
coefficientGrid = reshape(1:numCoefficients, basisSize);

lowerSubscripts = repmat({':'}, 1, numel(basisSize));
upperSubscripts = lowerSubscripts;
lowerSubscripts{dim} = 1:(basisSize(dim)-1);
upperSubscripts{dim} = 2:basisSize(dim);

lowerIndex = coefficientGrid(lowerSubscripts{:});
upperIndex = coefficientGrid(upperSubscripts{:});
lowerIndex = lowerIndex(:);
upperIndex = upperIndex(:);
numRows = numel(lowerIndex);

switch string(direction)
    case "increasing"
        rowValues = [ones(numRows,1); -ones(numRows,1)];
    case "decreasing"
        rowValues = [-ones(numRows,1); ones(numRows,1)];
    otherwise
        error('ConstrainedSpline:InvalidMonotonicDirection',  'direction must be "increasing" or "decreasing".');
end

rowIndex = [(1:numRows)'; (1:numRows)'];
columnIndex = [lowerIndex; upperIndex];
A = sparse(rowIndex, columnIndex, rowValues, numRows, numCoefficients);
