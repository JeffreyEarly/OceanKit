function Aeq = mesoscaleConstraintMatrix(psiKnotPoints, psiS, G, mesoscaleConstraint)
numReducedCoefficients = size(G, 2);
if mesoscaleConstraint == "none" || numReducedCoefficients == 0
    Aeq = sparse(0, numReducedCoefficients);
    return;
end

pointMatrix = constraintPointMatrix(psiKnotPoints, psiS);
Bqq = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[2 0 0]);
Brr = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 2 0]);

switch mesoscaleConstraint
    case "zeroVorticity"
        AeqFull = Bqq + Brr;
    case "zeroStrain"
        Bqr = TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 1 0]);
        AeqFull = [Bqr; Bqq - Brr];
    otherwise
        error("GriddedStreamfunction:UnknownMesoscaleConstraint", ...
            "Unknown mesoscale constraint '%s'.", mesoscaleConstraint);
end

Aeq = independentRows(AeqFull * G);
end

function pointMatrix = constraintPointMatrix(psiKnotPoints, psiS)
qSupport = supportPointsForDimension(psiKnotPoints{1}, psiS(1));
rSupport = supportPointsForDimension(psiKnotPoints{2}, psiS(2));
tSupport = supportPointsForDimension(psiKnotPoints{3}, psiS(3));
[qGrid, rGrid, tGrid] = ndgrid(qSupport, rSupport, tSupport);
pointMatrix = [qGrid(:), rGrid(:), tGrid(:)];
end

function supportPoints = supportPointsForDimension(knotVector, splineDegree)
uniqueKnots = unique(knotVector, "sorted");
intervalStarts = uniqueKnots(1:(end - 1));
intervalWidths = diff(uniqueKnots);
numInteriorPoints = splineDegree + 1;
interiorFractions = (1:numInteriorPoints).' / (numInteriorPoints + 1);
supportPoints = zeros(numInteriorPoints * nnz(intervalWidths > 0), 1);

nextIndex = 1;
for iInterval = 1:numel(intervalStarts)
    if intervalWidths(iInterval) <= 0
        continue;
    end

    intervalPoints = intervalStarts(iInterval) + interiorFractions * intervalWidths(iInterval);
    supportPoints(nextIndex:(nextIndex + numInteriorPoints - 1)) = intervalPoints;
    nextIndex = nextIndex + numInteriorPoints;
end
end

function Aeq = independentRows(Aeq)
numReducedCoefficients = size(Aeq, 2);
if isempty(Aeq) || numReducedCoefficients == 0
    Aeq = sparse(0, numReducedCoefficients);
    return;
end

[~, R, permutation] = qr(full(Aeq.'), "vector");
diagonal = abs(diag(R));
if isempty(diagonal)
    rankAeq = 0;
else
    tolerance = max(size(R)) * eps(class(R)) * max(diagonal);
    rankAeq = nnz(diagonal > tolerance);
end

if rankAeq == 0
    Aeq = sparse(0, numReducedCoefficients);
else
    rowIndices = sort(permutation(1:rankAeq));
    Aeq = sparse(Aeq(rowIndices, :));
end
end
