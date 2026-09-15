function self = createRecipe(options)
% Snapshot a family-owned continuous pairing descriptor.
arguments (Input)
    options.isAvailable (1,1) logical = true
    options.reason (1,1) string = ""
    options.variable (1,1) string = "u"
    options.columnLabels (1,:) string = strings(1,0)
    options.normalization (1,1) string = ""
    options.targetGramMatrix (:,:) double {mustBeFinite} = zeros(0)
    options.majorantGramMatrix (:,:) double {mustBeFinite} = zeros(0)
    options.supportsLeakage (1,1) logical = false
    options.supportsQuadratic (1,1) logical = false
    options.provenance (1,1) struct = struct()
    options.zDomain (1,2) double {mustBeReal,mustBeFinite} = [0 1]
    options.evaluateFunction (1,1) function_handle = @(z) zeros(numel(z),0)
    options.weightFunction (1,1) function_handle = @(z) zeros(size(z))
    options.spec (1,1) struct = struct()
end
arguments (Output)
    self (1,1) struct
end
self = options;
n = numel(self.columnLabels);
if self.isAvailable && (~isequal(size(self.targetGramMatrix),[n n]) || ~isequal(size(self.majorantGramMatrix),[n n]))
    error("IMProjection:InvalidTargetShape","Both target matrices must have one row and column per scientific label.");
end
if self.zDomain(1) >= self.zDomain(2)
    error("IMProjection:InvalidDomain","zDomain must be increasing.");
end
end
