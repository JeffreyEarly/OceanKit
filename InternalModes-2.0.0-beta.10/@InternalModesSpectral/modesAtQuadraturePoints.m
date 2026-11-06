function [F,G,h,zq] = modesAtQuadraturePoints(self,options)
% Return resolved modes on mode-adapted quadrature points.
%
% This helper chooses quadrature points from the roots of a resolved
% higher mode, temporarily evaluates the requested modes on those points,
% and restores the original output grid before returning.
%
% - Topic: Compute modes
% - Declaration: [F,G,h,zq] = modesAtQuadraturePoints(self,options)
% - Parameter self: InternalModesSpectral instance
% - Parameter options.nPoints: number of quadrature points and modes requested
% - Parameter options.omega: fixed frequency in radians per second; specify either `omega` or `k`
% - Parameter options.k: fixed horizontal wavenumber; specify either `omega` or `k`
% - Returns F: horizontal-velocity mode matrix on `zq`
% - Returns G: vertical-velocity mode matrix on `zq`
% - Returns h: equivalent-depth row vector
% - Returns zq: quadrature-point depth locations
arguments
    self InternalModesSpectral
    options.nPoints (1,1) double {mustBeInteger,mustBePositive} = 64
    options.omega double = []
    options.k double = []
end

hasOmega = ~isempty(options.omega);
hasK = ~isempty(options.k);
if hasOmega == hasK
    error('InternalModesSpectral:InvalidArguments', 'Specify exactly one of omega or k.');
end
if hasOmega && ~isscalar(options.omega)
    error('InternalModesSpectral:InvalidArguments', 'omega must be a scalar.');
end
if hasK && ~isscalar(options.k)
    error('InternalModesSpectral:InvalidArguments', 'k must be a scalar.');
end

originalZ = self.z;
cleanup = onCleanup(@() restoreOutputGrid(self, originalZ));

requestedModes = options.nPoints;
self.nEVP = max(self.nEVP, ceil(2.1*requestedModes));

maxIterations = 20;
for iIteration = 1:maxIterations
    if hasOmega
        self.gridFrequency = options.omega;
        [A,B] = self.eigenmatricesForFrequency(options.omega);
    else
        self.gridFrequency = 0;
        [A,B] = self.eigenmatricesForWavenumber(options.k);
    end

    [V_cheb,h] = self.solveGEP(A, B);
    resolvedModes = ceil(find(h > 0, 1, 'last')/2);
    if ~isempty(resolvedModes) && resolvedModes >= requestedModes
        zq = self.quadraturePointsForModes(requestedModes, V_cheb, h);
        self.z = zq;
        [F,G,h] = self.transformModesToSpatialDomain(V_cheb, h, requestedModes);
        return;
    end

    if isempty(resolvedModes) || resolvedModes == 0
        self.nEVP = ceil(1.5*self.nEVP);
    else
        self.nEVP = ceil(self.nEVP*(requestedModes/resolvedModes));
    end
    if self.shouldShowDiagnostics
        fprintf('Increasing nEVP to %d.\n', self.nEVP);
    end
end

error('InternalModesSpectral:NeedMorePoints', ...
    'Could not resolve %d modes after increasing nEVP to %d.', requestedModes, self.nEVP);
end

function restoreOutputGrid(self, originalZ)
self.z = originalZ;
end
