function errors = leakage(self,checkValues,checkNormSquared)
% Measure projected check profiles relative to their positive norms.
%
% Positive norms are supplied explicitly; no check modes are solved here.
% A zero-norm profile reports zero only when its sampled profile is zero;
% otherwise it reports infinity. Values are returned independently per input.
%
% - Topic: Measure projection quality
% - Declaration: errors = leakage(self,checkValues,checkNormSquared)
% - Parameter checkValues: finite real or complex sample-by-check matrix
% - Parameter checkNormSquared: positive squared reference norms, or exact zeros
% - Returns errors: row of relative projected norms
arguments (Input)
    self (1,1) IMProjection
    checkValues (:,:) double {mustBeFinite}
    checkNormSquared (:,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
arguments (Output)
    errors (1,:) double
end
if ~isvector(checkNormSquared) && ~isempty(checkNormSquared)
    error("IMProjection:InvalidShape","Squared norms must be a vector with one entry per input.");
end
checkNormSquared = reshape(checkNormSquared,1,[]);
if size(checkValues,2) ~= numel(checkNormSquared)
    error("IMProjection:InvalidShape","Supply one squared reference norm per check profile.");
end
coefficients = self.project(checkValues);
errors = IMProjection.relativeCoefficientNorm(coefficients,self.majorantGramMatrix,checkNormSquared);
zero = checkNormSquared == 0;
errors(zero) = Inf;
errors(zero & all(checkValues == 0,1)) = 0;
end
