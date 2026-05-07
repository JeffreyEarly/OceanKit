function z_g = quadraturePointsForModes(self,nPoints,G_cheb,h)
% Return quadrature points inferred from already-solved modes.
%
% - Topic: Developer topics
% - Developer: true
% - Declaration: z_g = quadraturePointsForModes(self,nPoints,G_cheb,h)
% - Parameter self: InternalModesSpectral instance
% - Parameter nPoints: number of quadrature points requested
% - Parameter G_cheb: sorted `G`-mode eigenvectors in Chebyshev coefficient space
% - Parameter h: sorted equivalent-depth vector
% - Returns z_g: depth locations of the quadrature points

resolvedModes = ceil(find(h > 0, 1, 'last')/2);
if isempty(resolvedModes) || resolvedModes < nPoints
    resolvedModeCount = 0;
    if ~isempty(resolvedModes)
        resolvedModeCount = resolvedModes;
    end
    error('InternalModesSpectral:NeedMorePoints', ...
        'Returned %d valid modes (%d quadrature points requested) using nEVP=%d.', ...
        resolvedModeCount, nPoints, self.nEVP);
end

if self.upperBoundary == UpperBoundary.mda
    rootMode = nPoints - 2;
elseif self.upperBoundary == UpperBoundary.rigidLid
    rootMode = nPoints - 1;
elseif self.upperBoundary == UpperBoundary.freeSurface && self.lowerBoundary == LowerBoundary.noSlip
    rootMode = nPoints - 1;
elseif self.upperBoundary == UpperBoundary.freeSurface
    rootMode = nPoints;
else
    rootMode = nPoints - 1;
end
rootMode = max(1, rootMode);

rootsVar = InternalModesSpectral.FindRootsFromChebyshevVector(G_cheb(:,rootMode), self.xDomain);
rootsVar(rootsVar < self.xMin) = self.xMin;
rootsVar(rootsVar > self.xMax) = self.xMax;
rootsVar = cat(1, self.xMin, rootsVar, self.xMax);
rootsVar = unique(rootsVar, 'stable');

while length(rootsVar) > nPoints
    rootsVar = sort(rootsVar);
    F = InternalModesSpectral.IntegrateChebyshevVector(G_cheb(:,rootMode));
    value = InternalModesSpectral.ValueOfFunctionAtPointOnGrid(rootsVar, self.xDomain, F);
    dv = diff(value);
    [~,minIndex] = min(abs(dv));
    rootsVar(minIndex+1) = [];
end

if length(rootsVar) < nPoints
    error('InternalModesSpectral:NeedMorePoints', ...
        'Returned %d unique roots (requested %d). Maybe need more EVP points.', length(rootsVar), nPoints);
end

z_g = reshape(rootsVar, [], 1);
z_g = InternalModesSpectral.fInverseBisection(self.x_function, z_g, min(self.zDomain), max(self.zDomain), 1e-12);
end
