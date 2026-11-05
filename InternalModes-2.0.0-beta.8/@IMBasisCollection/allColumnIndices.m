function columns = allColumnIndices(self,page)
% Resolve the default retained columns after validating the requested page.
arguments (Input)
    self (1,1) IMBasisCollection
    page (1,1) double {mustBeInteger,mustBePositive}
end
arguments (Output)
    columns (1,:) double
end
if page > numel(self.basisIndex)
    error("IMBasisCollection:InvalidPage","page must select an existing requested page.");
end
columns = 1:numel(self.metadata(self.basisIndex(page)).columnLabels);
end
