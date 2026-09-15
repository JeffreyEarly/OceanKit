function coefficients = project(self,values)
% Project sampled profiles with the signed physical pairing.
%
% - Topic: Apply projections
% - Declaration: coefficients = project(self,values)
% - Parameter values: finite real or complex sample-by-profile matrix
% - Returns coefficients: column-count-by-profile modal coefficients
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
coefficients = self.forwardMatrix*values;
end
