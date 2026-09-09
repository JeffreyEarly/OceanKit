function errors = productError(self,sampledProducts,referenceCoefficients,productNormSquared)
% Measure product projection errors with a positive coefficient norm.
%
% Reference coefficients must use the same signed target pairing, column
% coordinates, and normalization. A zero product norm returns zero only when
% sampled values and reference coefficients are both zero, and infinity
% otherwise. This method does not assert reference convergence or coverage.
%
% - Topic: Measure projection quality
% - Declaration: errors = productError(self,sampledProducts,referenceCoefficients,productNormSquared)
% - Parameter sampledProducts: finite real or complex sample-by-product values
% - Parameter referenceCoefficients: finite real or complex column-by-product coefficients
% - Parameter productNormSquared: nonnegative continuous squared product norms
% - Returns errors: row of relative coefficient-error magnitudes
arguments (Input)
    self (1,1) IMProjection
    sampledProducts (:,:) double {mustBeFinite}
    referenceCoefficients (:,:) double {mustBeFinite}
    productNormSquared (:,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
arguments (Output)
    errors (1,:) double
end
if ~isvector(productNormSquared) && ~isempty(productNormSquared)
    error("IMProjection:InvalidShape","Squared norms must be a vector with one entry per input.");
end
productNormSquared = reshape(productNormSquared,1,[]);
nProducts = size(sampledProducts,2);
if ~isequal(size(referenceCoefficients),[self.columnCount nProducts]) || numel(productNormSquared) ~= nProducts
    error("IMProjection:InvalidShape","Reference coefficients and norms must align with every product and projection column.");
end
if any(referenceCoefficients(~self.activeColumnMask,:) ~= 0,"all")
    error("IMProjection:InactiveReferenceCoefficient","Reference coefficients must be zero on inactive columns.");
end
difference = self.project(sampledProducts)-referenceCoefficients;
errors = IMProjection.relativeCoefficientNorm(difference,self.majorantGramMatrix,productNormSquared);
zero = productNormSquared == 0;
errors(zero) = Inf;
errors(zero & all(sampledProducts == 0,1) & all(referenceCoefficients == 0,1)) = 0;
end
