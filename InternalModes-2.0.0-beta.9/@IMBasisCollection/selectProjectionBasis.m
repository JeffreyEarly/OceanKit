function basis = selectProjectionBasis(self,page,columns)
% Resolve the requested basis and distinct columns for projection or assessment.
arguments (Input)
    self (1,1) IMBasisCollection
    page (1,1) double {mustBeInteger,mustBePositive}
    columns
end
if page > numel(self.basisIndex)
    error("IMBasisCollection:InvalidPage","page must select an existing requested page.");
end
if isequal(columns,":") || isequal(columns,':')
    columns = self.allColumnIndices(page);
end
if ~isa(columns,"double") || ~isrow(columns) && ~isempty(columns) || ~isreal(columns) ...
        || any(~isfinite(columns) | columns < 1 | columns ~= fix(columns),"all")
    error("IMProjection:InvalidColumns","columns must be a row of positive integer indices, [] or ':' for all columns.");
end
columns = reshape(columns,1,[]);
basis = self.bases{self.basisIndex(page)}.selectColumns(columns);
end
