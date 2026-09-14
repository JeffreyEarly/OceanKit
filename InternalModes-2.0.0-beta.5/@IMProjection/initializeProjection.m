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
    if self.projectionKind == "galerkin" && any(abs(targetNorms) <= 100*eps(max(1,abs(targetNorms))))
        error("IMProjection:InvalidTargetNorm","Every active target diagonal must be nonzero at numerical precision for Gram normalization.");
    end
    [~,failure] = chol(majorant(active,active));
    if failure ~= 0
        error("IMProjection:InvalidMajorant","majorantGramMatrix must be positive definite on active columns.");
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
self.forwardMatrix = zeros(nColumns,nSamples);
if isempty(active)
    self.sampledGramRank = 0;
    self.gramConditionNumber = NaN;
    self.targetGramIsPositiveDefinite = true;
    return
end
gramActive = sampleGram(active,active);
rankTolerance = max(size(gramActive))*eps(max(1,norm(gramActive,2)));
self.sampledGramRank = sum(svd(gramActive) > rankTolerance);
if self.sampledGramRank < numel(active)
    self.forwardMatrix(active,:) = pinv(gramActive,rankTolerance)*samplePairingMatrix(active,:);
else
    self.forwardMatrix(active,:) = gramActive \ samplePairingMatrix(active,:);
end
self.gramConditionNumber = cond(gramActive);
targetActive = targetGram(active,active);
self.targetGramIsPositiveDefinite = min(eig(targetActive)) > 100*eps(max(1,norm(targetActive,2)));
end

function matrix = symmetricMatrix(matrix,name)
if norm(matrix-matrix.',2) > 100*eps(max(1,norm(matrix,2)))
    error("IMProjection:NonSymmetricMatrix","%s must be symmetric.",name);
end
matrix = 0.5*(matrix+matrix.');
end
