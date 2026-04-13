function tKnot = normalizeKnotCell(tKnot, numDimensions)
% Normalize and validate a cell array of knot vectors.
if ~iscell(tKnot) || numel(tKnot) ~= numDimensions
    error('TensorSpline:InvalidKnotCell', 'knotPoints must be a cell array with one knot vector per dimension.');
end

tKnot = reshape(tKnot, 1, []);
for iDim = 1:numDimensions
    validateattributes(tKnot{iDim}, {'numeric'}, {'column','real','finite'});
    tKnot{iDim} = reshape(tKnot{iDim}, [], 1);
end
