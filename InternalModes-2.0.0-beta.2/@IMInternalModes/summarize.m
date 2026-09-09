function summarize(self, solver)
% Print a readable internal-mode EVP summary.
%
% `summarize` first prints the canonical EVP report, then appends the
% internal-mode interpretation: solved formulation, physical constants,
% equivalent-depth mapping, diagnostic-variable relation, factory-specific
% parameter names, and the available `F` and `G` inner products.
%
% ```matlab
% evp.summarize()
% ```
%
% prints output like:
%
% ```text
% Internal-mode context
%   formulation: G
%   f0: 0 s^-1
%   g: 9.81 m s^-2
%   equivalent depth: h = hFromEigenvalue(lambda)
%   mode family: hydrostatic
%   factory parameters: none
%
% Internal-mode variables
%   solved formulation: G
%   diagnostic variable: F
%   diagnostic relation: F_j(z) = h_j dG_j/dz(z)
%
% Internal-mode inner products
%   F
%     interior weight: 1
%     endpoint terms: none
%     inner product: known
%     reason: The hydrostatic value-only catalog gives an interior-only diagnostic inner product for this boundary combination.
%   G
%     interior weight: N^2(z)/g
%     endpoint terms: none
%     inner product: known
%     reason: The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.
%
% Internal-mode normalization
%   geostrophic normalization: available
%     shared scale: one factor per coupled (F,G) mode
%     convention: <G_j,G_j>_G = 1
%     implied: <F_j,F_j>_F = h_j
% ```
%
% For a fixed-frequency wave EVP, the diagnostic `F` inner product is
% reported as unavailable until a wave diagnostic catalog entry is added:
%
% ```text
% Internal-mode inner products
%   F
%     interior weight: 1
%     endpoint terms: none
%     inner product: unavailable
%     reason: The diagnostic inner product is available only for modeFamily="hydrostatic" EVPs in the value-only catalog.
% ```
%
% Mean-density-anomaly summaries instead report the surface-referenced
% `N2*G` integral for `F` and the signed generalized-energy normalization.
%
% - Topic: Summarize internal-mode EVPs
% - Declaration: summarize(evp,solver)
% - Parameter solver: optional solver for grid-level coefficient and mode-selection assessment
arguments
    self IMInternalModes
    solver = []
end

if isempty(solver)
    summarize@IMEigenvalueProblem(self);
else
    summarize@IMEigenvalueProblem(self, solver);
end

fprintf('\nInternal-mode context\n');
fprintf('  formulation: %s\n', self.formulation);
fprintf('  f0: %s s^-1\n', formatNumber(self.f0));
fprintf('  g: %s m s^-2\n', formatNumber(self.g));
fprintf('  equivalent depth: h = hFromEigenvalue(lambda)\n');
fprintf('  mode family: %s\n', self.modeFamily);
fprintf('  factory parameters: %s\n', factoryParameterText(self.parameters));

fprintf('\nInternal-mode variables\n');
fprintf('  solved formulation: %s\n', self.formulation);
fprintf('  diagnostic variable: %s\n', diagnosticVariable(self));
fprintf('  diagnostic relation: %s\n', diagnosticRelationText(self));

fprintf('\nInternal-mode inner products\n');
printInnerProductSpec(self.innerProduct("F"));
printInnerProductSpec(self.innerProduct("G"));

fprintf('\nInternal-mode normalization\n');
printNormalizationSummary(self);
end

function text = factoryParameterText(parameters)
names = setdiff(string(fieldnames(parameters)), ["f0", "g", "formulation"], "stable");
if isempty(names)
    text = "none";
else
    text = join(names, ", ");
end
end

function variable = diagnosticVariable(evp)
if evp.formulation == "G"
    variable = "F";
else
    variable = "G";
end
end

function text = diagnosticRelationText(evp)
if evp.modeFamily == "meanDensityAnomaly"
    text = "F_j(z) = (1/g) integral_z^surface N^2(z') G_j(z') dz'";
elseif evp.formulation == "G"
    text = "F_j(z) = h_j dG_j/dz(z)";
elseif evp.modeFamily == "hydrostatic"
    text = "G_j(z) = -g/N^2(z) dF_j/dz(z)";
else
    text = "G_j(z) = GfromFz(z,dF_j/dz,h_j,ctx)";
end
end

function printNormalizationSummary(evp)
if evp.modeFamily == "hydrostatic"
    fprintf('  geostrophic normalization: available\n');
    fprintf('    shared scale: one factor per coupled (F,G) mode\n');
    fprintf('    convention: <G_j,G_j>_G = 1\n');
    fprintf('    implied: <F_j,F_j>_F = h_j\n');
elseif evp.modeFamily == "meanDensityAnomaly"
    fprintf('  signed generalized-energy normalization: available\n');
    fprintf('    shared scale: one positive factor per aligned (F,G) mode\n');
    fprintf('    convention: <G_j,G_j>_G = epsilon_j, epsilon_j in {-1,+1}\n');
else
    fprintf('  geostrophic normalization: unavailable for this mode family\n');
end
end

function printInnerProductSpec(spec)
fprintf('  %s\n', spec.variable);
fprintf('    interior weight: %s\n', interiorWeightText(spec.variable));
fprintf('    endpoint terms: %s\n', endpointTermsText(spec));
if spec.hasInnerProduct
    fprintf('    inner product: known\n');
else
    fprintf('    inner product: unavailable\n');
end
fprintf('    reason: %s\n', spec.reason);
end

function text = interiorWeightText(variable)
if string(variable) == "F"
    text = "1";
else
    text = "N^2(z)/g";
end
end

function text = endpointTermsText(spec)
weights = [spec.surfaceWeights; spec.bottomWeights];
terms = spec.endpointInnerProductTerms;
if isempty(weights)
    weightParts = strings(1,0);
else
    weightParts = strings(1,numel(weights));
    for iWeight = 1:numel(weights)
        weight = weights(iWeight);
        expression = endpointExpression(weight.c, -weight.d, weight.location);
        weightParts(iWeight) = weight.location + ": " + formatSignedNumber(weight.coefficient) + " * (" + expression + ")^2";
    end
end

if isempty(terms)
    termParts = strings(1,0);
else
    termParts = strings(1,numel(terms));
    for iTerm = 1:numel(terms)
        term = terms(iTerm);
        termParts(iTerm) = term.location + ": " ...
            + formatSignedNumber(term.coefficient) + " * " ...
            + term.variable + "(" + term.location + ")^2";
    end
end

parts = [weightParts termParts];
if isempty(parts)
    text = "none";
else
    text = join(parts, "; ");
end
end

function expression = endpointExpression(valueCoefficient, fluxCoefficient, location)
labels = ["u(" + location + ")", "p(" + location + ")*du/dz(" + location + ")"];
expression = signedExpression([valueCoefficient fluxCoefficient], labels);
end

function expression = signedExpression(coefficients, labels)
expression = "";
tolerance = 100*eps*max(1,max(abs(coefficients)));
for iCoefficient = 1:numel(coefficients)
    coefficient = coefficients(iCoefficient);
    if abs(coefficient) <= tolerance
        continue;
    end
    term = formatTerm(abs(coefficient), labels(iCoefficient));
    if strlength(expression) == 0
        if coefficient < 0
            expression = "-" + term;
        else
            expression = term;
        end
    elseif coefficient < 0
        expression = expression + " - " + term;
    else
        expression = expression + " + " + term;
    end
end
if strlength(expression) == 0
    expression = "0";
end
end

function text = formatTerm(coefficient, label)
if coefficient == 1
    text = label;
else
    text = formatNumber(coefficient) + "*" + label;
end
end

function text = formatNumber(value)
if value == 0
    text = "0";
else
    text = string(sprintf('%.6g', value));
end
end

function text = formatSignedNumber(value)
if value == 0
    text = "+0";
else
    text = string(sprintf('%+.6g', value));
end
end
