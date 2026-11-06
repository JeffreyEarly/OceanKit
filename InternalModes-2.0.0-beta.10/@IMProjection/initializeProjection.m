function self = initializeProjection(self,samplePairingMatrix,sampleGram,targetGram,majorant,activeMask,labels)
% Initialize shared coefficient-space state from validated finite arrays.
nColumns = size(samplePairingMatrix,1);
nSamples = size(samplePairingMatrix,2);
if ~isequal(size(sampleGram),[nColumns nColumns]) || ~isequal(size(targetGram),[nColumns nColumns]) || numel(activeMask) ~= nColumns || numel(labels) ~= nColumns
    error("IMProjection:InvalidShape","Coefficient systems, active mask, and labels must match the pairing row count.");
end
if any(ismissing(labels) | strlength(labels) == 0) || numel(unique(labels)) ~= nColumns
    error("IMProjection:InvalidLabels","Column labels must be nonempty, nonmissing, and distinct.");
end
sampleGram = symmetricMatrix(sampleGram,"sampleGram");
targetGram = symmetricMatrix(targetGram,"targetGramMatrix");
if isempty(majorant)
    if ~isequal(targetGram,diag(diag(targetGram)))
        error("IMProjection:MajorantRequired","A nondiagonal target Gram matrix requires an explicit positive majorantGramMatrix.");
    end
    majorant = diag(abs(diag(targetGram)));
end
if ~isequal(size(majorant),[nColumns nColumns])
    error("IMProjection:InvalidShape","majorantGramMatrix must have one row and column per output coordinate.");
end
majorant = symmetricMatrix(majorant,"majorantGramMatrix");
active = find(activeMask);
for iColumn = find(~activeMask)
    values = [samplePairingMatrix(iColumn,:).';sampleGram(:,iColumn);targetGram(:,iColumn);majorant(:,iColumn)];
    if norm(values,2) > 1e3*eps(max(1,norm(values,2)))
        error("IMProjection:InvalidInactiveColumn","Inactive coordinates must have numerically zero pairings and coefficient metrics.");
    end
end
if ~isempty(active)
    targetNorms = diag(targetGram(active,active));
    if self.projectionKind == "galerkin" && any(targetNorms == 0)
        error("IMProjection:InvalidTargetNorm","Every active target diagonal must be nonzero relative to its positive coefficient metric for Gram normalization.");
    end
    majorantDiagonal = diag(majorant(active,active));
    if any(majorantDiagonal <= 0)
        error("IMProjection:InvalidMajorant","majorantGramMatrix must be positive definite on active columns.");
    end
    coefficientScale = sqrt(majorantDiagonal);
    normalizedMajorant = (majorant(active,active)./coefficientScale)./coefficientScale.';
    [majorantFactor,failure] = chol(normalizedMajorant);
    if failure ~= 0
        error("IMProjection:InvalidMajorant","majorantGramMatrix must be positive definite on active columns.");
    end
    relativeTargetNorms = abs(targetNorms)./majorantDiagonal;
    if self.projectionKind == "galerkin" && any(relativeTargetNorms <= 100*eps(max(1,relativeTargetNorms)))
        error("IMProjection:InvalidTargetNorm","Every active target diagonal must be nonzero relative to its positive coefficient metric for Gram normalization.");
    end
end
self.sampleCount = nSamples;
self.columnCount = nColumns;
self.samplePairingMatrix = samplePairingMatrix;
self.gramMatrix = sampleGram;
self.targetGramMatrix = targetGram;
self.majorantGramMatrix = majorant;
self.activeColumnMask = activeMask;
self.columnLabels = labels;
self.recovery = [];
self.physicalDualMatrix = [];
self.physicalDualStatus = "unavailable";
if isempty(active)
    self.coefficientScale_ = zeros(0);
    self.physicalDualMatrix = zeros(nColumns,nSamples);
    self.physicalDualStatus = "available";
    self.sampledGramRank = 0;
    self.gramConditionNumber = NaN;
    self.targetGramIsPositiveDefinite = true;
    return
end
gramActive = sampleGram(active,active);
% Factor and test the signed pairing in dimensionless coefficient coordinates.
% The majorant supplies units, not a replacement physical pairing. In these
% coordinates column rescaling cannot turn a weak but independent mode into
% a numerical null column.
normalizedGram = (gramActive./coefficientScale)./coefficientScale.';
normalizedGram = majorantFactor.'\normalizedGram/majorantFactor;
normalizedGram = 0.5*(normalizedGram+normalizedGram.');
normalizedPairing = majorantFactor.'\(samplePairingMatrix(active,:)./coefficientScale);
self.coefficientScale_ = diag(1./coefficientScale)/majorantFactor;
rankTolerance = max(size(normalizedGram))*eps(norm(normalizedGram,2));
self.sampledGramRank = sum(svd(normalizedGram) > rankTolerance);
self.gramConditionNumber = cond(gramActive);
targetActive = targetGram(active,active);
normalizedTarget = (targetActive./coefficientScale)./coefficientScale.';
normalizedTarget = majorantFactor.'\normalizedTarget/majorantFactor;
normalizedTarget = 0.5*(normalizedTarget+normalizedTarget.');
self.targetGramIsPositiveDefinite = min(eig(normalizedTarget)) > 100*eps(norm(normalizedTarget,2));
if self.projectionKind == "prescribedDual"
    % This is the caller's physical coefficient normalization, not a sampled
    % synthesis inverse. A singular supplied system is explicitly unavailable.
    if self.sampledGramRank < numel(active)
        self.physicalDualStatus = "singularCoefficientSystem";
    else
        self.physicalDualMatrix = zeros(nColumns,nSamples);
        self.physicalDualMatrix(active,:) = self.coefficientScale_*(normalizedGram\normalizedPairing);
        if any(~isfinite(self.physicalDualMatrix),"all") || any(any(samplePairingMatrix(active,:) ~= 0,1) & all(self.physicalDualMatrix(active,:) == 0,1))
            self.physicalDualMatrix = [];
            self.physicalDualStatus = "unrepresentableCoefficientSystem";
        else
            self.physicalDualStatus = "available";
        end
    end
else
    targetTolerance = max(size(normalizedTarget))*eps(norm(normalizedTarget,2));
    if sum(svd(normalizedTarget) > targetTolerance) < numel(active)
        error("IMProjection:SingularTarget","A physical coefficient dual requires a nonsingular continuous target on the active columns.");
    end
    self.physicalDualMatrix = zeros(nColumns,nSamples);
    self.physicalDualMatrix(active,:) = self.coefficientScale_*(normalizedTarget\normalizedPairing);
    if any(~isfinite(self.physicalDualMatrix),"all") || any(any(samplePairingMatrix(active,:) ~= 0,1) & all(self.physicalDualMatrix(active,:) == 0,1))
        self.physicalDualMatrix = [];
        self.physicalDualStatus = "unrepresentableCoefficientSystem";
    else
        self.physicalDualStatus = "available";
    end
end
end

function matrix = symmetricMatrix(matrix,name)
if norm(matrix-matrix.',2) > 100*eps(max(1,norm(matrix,2)))
    error("IMProjection:NonSymmetricMatrix","%s must be symmetric.",name);
end
matrix = 0.5*(matrix+matrix.');
end
