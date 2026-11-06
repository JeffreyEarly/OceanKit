function errors = relativeCoefficientNorm(coefficients,majorantGramMatrix,normSquared)
% Measure coefficient magnitudes in an already validated positive metric.
%
% Internal assessment kernels supply their validated active majorant and
% reference norms. This helper owns the shared Hermitian quadratic form and
% ratio, not projection, reference qualification, or zero-profile detection.
% A zero denominator with zero numerator returns NaN until the owning caller
% classifies the source profile. Roundoff-sized negative squared norms clip
% to zero; materially negative values reject an invalid metric calculation.
%
% - Topic: Measure projection quality
% - Developer: true
% - Declaration: errors = IMProjection.relativeCoefficientNorm(coefficients,majorantGramMatrix,normSquared)
arguments (Input)
    coefficients (:,:) double {mustBeFinite}
    majorantGramMatrix (:,:) double {mustBeReal,mustBeFinite}
    normSquared (:,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
arguments (Output)
    errors (1,:) double
end
if ~isequal(size(majorantGramMatrix),[size(coefficients,1) size(coefficients,1)]) || numel(normSquared) ~= size(coefficients,2) || (~isvector(normSquared) && ~isempty(normSquared))
    error("IMProjection:InvalidShape","Coefficient metrics and reference norms must match the coefficient array.");
end
weighted = majorantGramMatrix*coefficients;
numeratorSquared = real(sum(conj(coefficients).*weighted,1));
roundoffScale = sum(abs(coefficients).*abs(weighted),1);
if any(numeratorSquared < -100*eps(max(1,roundoffScale)))
    error("IMProjection:InvalidMajorant","The supplied coefficient metric produced a materially negative squared norm.");
end
errors = sqrt(max(0,numeratorSquared)./reshape(normSquared,1,[]));
end
