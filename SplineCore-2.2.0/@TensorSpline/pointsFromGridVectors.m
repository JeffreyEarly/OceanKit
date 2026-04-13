function [pointMatrix, gridSize] = pointsFromGridVectors(gridVectors)
% Convert rectilinear grid vectors into an explicit point matrix.
%
% Use this helper to convert rectilinear grid vectors into the
% point-matrix format accepted by `TensorSpline.matrixForPointMatrix`.
%
% If `[X1,...,Xd] = ndgrid(gridVectors{:})`, then the returned matrix is
%
% $$
% \mathrm{pointMatrix} = [X_1(:)\ \cdots\ X_d(:)].
% $$
%
% ```matlab
% [points, gridSize] = TensorSpline.pointsFromGridVectors({x,y});
% B = TensorSpline.matrixForPointMatrix(points, knotPoints=tKnot, S=[3 3]);
% ```
%
% - Topic: Build spline bases
% - Declaration: [pointMatrix,gridSize] = pointsFromGridVectors(gridVectors)
% - Parameter gridVectors: cell array of grid vectors
% - Returns pointMatrix: matrix with one row per grid point
% - Returns gridSize: number of points along each dimension
arguments
    gridVectors cell
end

numDimensions = numel(gridVectors);
gridVectors = TensorSpline.normalizeKnotCell(gridVectors, numDimensions);
gridSize = cellfun(@numel, gridVectors);

grids = cell(1, numDimensions);
[grids{:}] = ndgrid(gridVectors{:});

pointMatrix = zeros(prod(gridSize), numDimensions);
for iDim = 1:numDimensions
    pointMatrix(:,iDim) = grids{iDim}(:);
end
