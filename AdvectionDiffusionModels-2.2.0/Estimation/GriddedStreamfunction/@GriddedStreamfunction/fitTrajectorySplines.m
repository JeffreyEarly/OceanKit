function fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS, mesoscaleConstraint, buildDecomposition, sampleData)
if isempty(sampleData)
    sampledData = GriddedStreamfunction.sampleTrajectoryData(trajectories);
else
    sampledData = sampleData;
end
tCell = sampledData.tCell;
allT = sampledData.allT;
allX = sampledData.allX;
allY = sampledData.allY;
allXDot = sampledData.allXDot;
allYDot = sampledData.allYDot;
fitSupportTimes = sampledData.sampleSupportTimes;

if isempty(fastKnotPoints)
    representativeTimes = GriddedStreamfunction.representativeObservationTimes(tCell);
    % Anchor the automatic fast basis to the full observation interval so
    % asynchronous deployment gaps do not clip the spline domain.
    fastSupportTimes = unique([fitSupportTimes(1); representativeTimes; fitSupportTimes(end)], "sorted");
    if numel(fastSupportTimes) < fastS + 1
        error("GriddedStreamfunction:InsufficientFastSupportTimes", ...
            "The representative times must contain at least fastS + 1 unique samples.");
    end

    fastKnotPoints = BSpline.knotPointsForDataPoints( ...
        fastSupportTimes, S=fastS, splineDOF=numel(fastSupportTimes));
else
    representativeTimes = zeros(0, 1);
    if any(allT < fastKnotPoints(1) | allT > fastKnotPoints(end))
        error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
            "Observation times must lie inside the supplied fast spline domain.");
    end
end

fastState = GriddedStreamfunction.fastBasisState(allT, fastKnotPoints, fastS, true);
centerCoefficients = fastState.solver \ [allX, allY];
bX = centerCoefficients(:, 1);
bY = centerCoefficients(:, 2);
centerXSpline = TensorSpline.fromKnotPoints(fastKnotPoints, bX, S=fastS);
centerYSpline = TensorSpline.fromKnotPoints(fastKnotPoints, bY, S=fastS);
centerOfMassTrajectory = TrajectorySpline(t=fitSupportTimes, x=centerXSpline, y=centerYSpline);

mxAll = fastState.B * bX;
myAll = fastState.B * bY;
% \dot{m} comes from differentiating the fitted fast COM spline.
mxDotAll = fastState.dB * bX;
myDotAll = fastState.dB * bY;
qAll = allX - mxAll;
rAll = allY - myAll;

if isempty(psiKnotPoints)
    psiKnotPoints = GriddedStreamfunction.defaultPsiKnotPoints(qAll, rAll, allT, psiS);
end

psiBasisSize = TensorSpline.basisSizeFromKnotCell(psiKnotPoints, psiS + 1);
qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
tDomain = [psiKnotPoints{3}(1), psiKnotPoints{3}(end)];
if any(qAll < qDomain(1) | qAll > qDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered q observations must lie inside the supplied psi spline domain.");
end
if any(rAll < rDomain(1) | rAll > rDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered r observations must lie inside the supplied psi spline domain.");
end
if any(allT < tDomain(1) | allT > tDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Observation times must lie inside the supplied psi time domain.");
end

pointMatrix = [qAll, rAll, allT];
Hx = -sparse(TensorSpline.matrixForPointMatrix( ...
    pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 1 0]));
Hy = sparse(TensorSpline.matrixForPointMatrix( ...
    pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 0 0]));

numSpatialCoefficients = prod(psiBasisSize(1:2));
if numSpatialCoefficients == 1
    Gspatial = sparse(1, 0);
else
    spatialSupportCounts = max(2 * psiBasisSize(1:2), psiBasisSize(1:2) + 2);
    qGaugePoints = linspace(qDomain(1), qDomain(2), spatialSupportCounts(1)).';
    rGaugePoints = linspace(rDomain(1), rDomain(2), spatialSupportCounts(2)).';
    [qGaugeGrid, rGaugeGrid] = ndgrid(qGaugePoints, rGaugePoints);
    spatialPointMatrix = [qGaugeGrid(:), rGaugeGrid(:)];
    spatialBasisMatrix = TensorSpline.matrixForPointMatrix( ...
        spatialPointMatrix, knotPoints=psiKnotPoints(1:2), S=psiS(1:2));
    constantGaugeVector = spatialBasisMatrix \ ones(size(qGaugeGrid(:)));
    constantGaugeResidual = spatialBasisMatrix * constantGaugeVector - ones(size(qGaugeGrid(:)));
    if norm(constantGaugeResidual, inf) > 1e-10
        error("GriddedStreamfunction:InvalidPsiGauge", ...
            "The spatial spline basis must represent a constant streamfunction gauge.");
    end

    Gspatial = sparse(null(reshape(constantGaugeVector, [], 1).', "r"));
end
G = kron(speye(psiBasisSize(3)), Gspatial);

Hx = Hx * G;
Hy = Hy * G;
% Reuse the fast-basis factorization in the same least-squares projector B * (B \ A).
McomX = sparse(fastState.B * (fastState.solver \ Hx));
McomY = sparse(fastState.B * (fastState.solver \ Hy));
Mcom = [McomX; McomY];
Mrel = [Hx - McomX; Hy - McomY];
dcom = [mxDotAll; myDotAll];
drel = [allXDot - mxDotAll; allYDot - myDotAll];
Aeq = GriddedStreamfunction.mesoscaleConstraintMatrix(psiKnotPoints, psiS, G, mesoscaleConstraint);
mesoscaleDegreesOfFreedom = size(G, 2) - size(Aeq, 1);
X = [Mcom; Mrel];
d = [dcom; drel];
% d = [Mcom; Mrel] * alpha.
if isempty(Aeq)
    alpha = X \ d;
else
    alpha = ConstrainedSpline.constrainedWeightedSolution(X.' * X, X.' * d, Aeq, zeros(size(Aeq, 1), 1), [], []);
end
streamfunctionCoefficients = G * alpha;

streamfunctionSpline = TensorSpline.fromKnotPoints(psiKnotPoints, streamfunctionCoefficients, S=psiS);

self.streamfunctionSpline = streamfunctionSpline;
self.observedTrajectories = trajectories;
self.observedTrajectorySampleData = sampledData;
self.centerOfMassTrajectory = centerOfMassTrajectory;
self.mesoscaleConstraint = mesoscaleConstraint;
self.representativeTimes = representativeTimes;
self.fitSupportTimes = fitSupportTimes;
self.mesoscaleDegreesOfFreedom = mesoscaleDegreesOfFreedom;
if buildDecomposition
    [backgroundTrajectory, decomposition] = decomposeTrajectorySet(self, trajectories);
    self.backgroundTrajectory = backgroundTrajectory;
    self.fixedFrameBackgroundTrajectories = decomposition.fixedFrame.background;
    self.fixedFrameMesoscaleTrajectories = decomposition.fixedFrame.mesoscale;
    self.fixedFrameSubmesoscaleTrajectories = decomposition.fixedFrame.submesoscale;
    self.centeredFrameMesoscaleTrajectories = decomposition.centeredFrame.mesoscale;
    self.centeredFrameSubmesoscaleTrajectories = decomposition.centeredFrame.submesoscale;
else
    self.backgroundTrajectory = TrajectorySpline.empty(0, 1);
    self.fixedFrameBackgroundTrajectories = TrajectorySpline.empty(0, 1);
    self.fixedFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1);
    self.fixedFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1);
    self.centeredFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1);
    self.centeredFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1);
end
end
