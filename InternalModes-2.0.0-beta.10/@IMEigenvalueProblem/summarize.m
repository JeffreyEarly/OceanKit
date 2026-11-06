function summarize(self, solver)
% Print a readable mathematical summary of this EVP.
%
% `summarize` prints the canonical scalar equation, physical domain,
% boundary conditions, endpoint norm weights, and parameter names. When a
% solver is supplied, it also samples the coefficients on the solver grid,
% assembles the left matrix for the zero-mode check, and prints the
% grid-level negative and zero mode assessment. It does not solve the EVP.
%
% ```matlab
% evp.summarize()
% ```
%
% prints output like:
%
% ```text
% dirichlet
%
% Canonical EVP
%   -(d/dz)(p(z) du/dz) + q(z) u = lambda r(z) u
%   z in [-1, 0]
%
% Boundary conditions
%   surface: u(surface) = 0
%   bottom: u(bottom) = 0
% ```
%
% ```matlab
% solver = IMSolverSpectral(nEVP=64);
% evp.summarize(solver)
% ```
%
% prints additional solver-grid output like:
%
% ```text
% Solver assessment
%   solver: IMSolverSpectral
%   coordinate: z
%   grid size: 64
%
% Coefficient ranges on solver grid
%   p: [1, 1]
%   q: [0, 0]
%   r: [1, 1]
%
% Mode selection assessment
%   negative modes: expected 0
%   zero mode: absent
% ```
%
% - Topic: Summarize EVPs
% - Declaration: summarize(evp,solver)
% - Parameter solver: optional solver for grid-level coefficient and mode-selection assessment
arguments
    self IMEigenvalueProblem
    solver = []
end

if ~isempty(solver) && ~isa(solver, "IMSolver")
    error("IMEigenvalueProblem:InvalidSolver", ...
        "The optional summarize input must be an IMSolver instance.");
end

fprintf('%s\n\n', self.name);
fprintf('Canonical EVP\n');
fprintf('  -(d/dz)(p(z) du/dz) + q(z) u = lambda r(z) u\n');
fprintf('  z in [%s, %s]\n\n', formatNumber(self.zDomain(1)), formatNumber(self.zDomain(2)));

fprintf('Boundary conditions\n');
printBoundary("surface", self.surfaceBoundary);
printBoundary("bottom", self.bottomBoundary);

fprintf('\nEndpoint norm weights\n');
printEndpointWeights(self.endpointWeights());

fprintf('\nParameters\n');
printParameterNames(self.parameters);

if ~isempty(solver)
    configuredSolver = solver.configuredForEVP(self);
    context = self.contextForSolver(configuredSolver);
    z = configuredSolver.zNative(:);
    [pValues, qValues, rValues] = self.coefficientValues(z, context);
    [A, ~] = self.assemble(configuredSolver);
    diagnostics = self.modeSelectionDiagnostics(configuredSolver, A);
    printSolverAssessment(configuredSolver, pValues, qValues, rValues, diagnostics);
end
end

function printBoundary(location, boundary)
location = string(location);
leftExpression = endpointExpression(boundary.a, -boundary.b, location);
rightExpression = endpointExpression(boundary.c, -boundary.d, location);

if ~boundary.isEigenvalueDependent()
    if boundary.b == 0
        fprintf('  %s: u(%s) = 0\n', location, location);
    elseif boundary.a == 0
        fprintf('  %s: p(%s)*du/dz(%s) = 0\n', location, location, location);
    else
        fprintf('  %s: %s = 0\n', location, leftExpression);
    end
else
    fprintf('  %s: -[%s] = lambda*[%s]\n', location, leftExpression, rightExpression);
end
end

function printEndpointWeights(weights)
if isempty(weights)
    fprintf('  none\n');
    return;
end

for iWeight = 1:numel(weights)
    weight = weights(iWeight);
    expression = endpointExpression(weight.c, -weight.d, weight.location);
    fprintf('  %s: %s * (%s)^2\n', weight.location, formatSignedNumber(weight.coefficient), expression);
end
end

function printParameterNames(parameters)
names = string(fieldnames(parameters));
if isempty(names)
    fprintf('  none\n');
else
    fprintf('  names: %s\n', join(names, ", "));
end
end

function printSolverAssessment(configuredSolver, pValues, qValues, rValues, diagnostics)
fprintf('\nSolver assessment\n');
fprintf('  solver: %s\n', class(configuredSolver));
if isprop(configuredSolver, "coordinateKind")
    fprintf('  coordinate: %s\n', configuredSolver.coordinateKind);
else
    fprintf('  coordinate: finiteDifference\n');
end
if isprop(configuredSolver, "nEVP")
    fprintf('  grid size: %d\n', configuredSolver.nEVP);
end

fprintf('\nCoefficient ranges on solver grid\n');
fprintf('  p: %s\n', formatRange(pValues));
fprintf('  q: %s\n', formatRange(qValues));
fprintf('  r: %s\n', formatRange(rValues));

fprintf('\nMode selection assessment\n');
fprintf('  negative modes: %s\n', negativeModeText(diagnostics));
fprintf('  zero mode: %s\n', zeroModeText(diagnostics));
fprintf('  reason: %s\n', diagnostics.reason);
end

function text = negativeModeText(diagnostics)
minCount = diagnostics.minNegativeEigenvalueCount;
maxCount = diagnostics.maxNegativeEigenvalueCount;
if ~isnumeric(maxCount)
    text = "unknown";
elseif minCount == maxCount
    text = "expected " + string(minCount);
else
    text = "possible range [" + string(minCount) + ", " + string(maxCount) + "]";
end
end

function text = zeroModeText(diagnostics)
text = diagnostics.zeroModeStatus;
if diagnostics.zeroModeStatus == "present"
    text = text + " (count " + string(diagnostics.zeroModeCount) + ")";
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

function text = formatRange(values)
values = values(:);
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    text = "all nonfinite";
    return;
end
text = "[" + formatNumber(min(finiteValues)) + ", " + formatNumber(max(finiteValues)) + "]";
nNonfinite = nnz(~isfinite(values));
if nNonfinite > 0
    text = text + " (" + string(nNonfinite) + " nonfinite)";
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
