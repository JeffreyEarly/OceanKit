function coefficients = projectPhysical(self,values)
% Apply the physical coefficient dual without sampled-overlap correction.
%
% For a Galerkin physical record, this is the continuous target normalization
% applied to the raw sampled pairings. Quadrature errors remain visible.
% A prescribed source uses its explicitly supplied coefficient system.
% Physical projection availability does not require successful state recovery.
%
% - Topic: Apply projections
% - Declaration: coefficients = projectPhysical(projection,values)
% - Parameter values: finite real or complex sample-by-profile values
% - Returns coefficients: physical coefficient estimates in the declared coordinates
arguments (Input)
    self (1,1) IMProjection
    values (:,:) double {mustBeFinite}
end
arguments (Output)
    coefficients (:,:) double
end
if size(values,1) ~= self.sampleCount
    error("IMProjection:InvalidShape","values must have one row per sample.");
end
if self.physicalDualStatus ~= "available"
    error("IMProjection:UnavailablePhysicalDual","The physical coefficient system is unavailable (%s); raw pairings remain available.",self.physicalDualStatus);
end
coefficients = self.physicalDualMatrix*values;
end
