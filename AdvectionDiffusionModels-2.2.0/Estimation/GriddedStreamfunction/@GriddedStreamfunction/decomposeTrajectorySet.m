function [backgroundTrajectory, decomposition] = decomposeTrajectorySet(self, trajectories)
arguments (Input)
    self (1,1) GriddedStreamfunction
    trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
end
arguments (Output)
    backgroundTrajectory (1,1) TrajectorySpline
    decomposition struct
end

trajectories = reshape(trajectories, [], 1);
sampleData = decompositionSampleData(self, trajectories);
nTrajectories = numel(sampleData.tCell);

backgroundXCoefficients = sampleData.B \ sampleData.uBackgroundObserved;
backgroundYCoefficients = sampleData.B \ sampleData.vBackgroundObserved;
backgroundXSpline = TensorSpline.fromKnotPoints(self.fastKnotPoints, backgroundXCoefficients, S=self.fastS);
backgroundYSpline = TensorSpline.fromKnotPoints(self.fastKnotPoints, backgroundYCoefficients, S=self.fastS);

backgroundTrajectoryXSpline = cumsum(backgroundXSpline);
backgroundTrajectoryYSpline = cumsum(backgroundYSpline);
backgroundTrajectoryXSpline = TensorSpline.fromKnotPoints( ...
    backgroundTrajectoryXSpline.knotPoints, backgroundTrajectoryXSpline.xi, ...
    S=backgroundTrajectoryXSpline.S, xMean=-backgroundTrajectoryXSpline(sampleData.sampleSupportTimes(1)));
backgroundTrajectoryYSpline = TensorSpline.fromKnotPoints( ...
    backgroundTrajectoryYSpline.knotPoints, backgroundTrajectoryYSpline.xi, ...
    S=backgroundTrajectoryYSpline.S, xMean=-backgroundTrajectoryYSpline(sampleData.sampleSupportTimes(1)));
backgroundTrajectory = TrajectorySpline( ...
    t=sampleData.sampleSupportTimes, x=backgroundTrajectoryXSpline, y=backgroundTrajectoryYSpline);

backgroundTrajectories = trajectories;
mesoscaleTrajectories = trajectories;
submesoscaleTrajectories = trajectories;
centeredMesoscaleTrajectories = trajectories;
centeredSubmesoscaleTrajectories = trajectories;
cachedTimes = cell(0, 1);
cachedDegrees = zeros(0, 1);
cachedBasisMatrices = cell(0, 1);
cachedKnotPoints = cell(0, 1);

sampleStartIndex = 1;
for iTrajectory = 1:nTrajectories
    ti = sampleData.tCell{iTrajectory};
    nSamples = numel(ti);
    sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
    componentS = min(3, nSamples - 1);
    x0 = sampleData.allX(sampleIndices(1));
    y0 = sampleData.allY(sampleIndices(1));
    q0 = sampleData.qAll(sampleIndices(1));
    r0 = sampleData.rAll(sampleIndices(1));

    backgroundTrajectories(iTrajectory) = reanchoredBackgroundTrajectory(backgroundTrajectory, ti);
    [mesoscaleTrajectories(iTrajectory), cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledVelocityTrajectory( ...
        ti, sampleData.uMesoscaleObserved(sampleIndices), sampleData.vMesoscaleObserved(sampleIndices), componentS, x0, y0, ...
        cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints);
    [submesoscaleTrajectories(iTrajectory), cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledVelocityTrajectory( ...
        ti, sampleData.uSubmesoscaleObserved(sampleIndices), sampleData.vSubmesoscaleObserved(sampleIndices), componentS, 0, 0, ...
        cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints);
    [centeredMesoscaleTrajectories(iTrajectory), cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledVelocityTrajectory( ...
        ti, sampleData.uMesoscaleRelativeObserved(sampleIndices), sampleData.vMesoscaleRelativeObserved(sampleIndices), componentS, q0, r0, ...
        cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints);
    [centeredSubmesoscaleTrajectories(iTrajectory), cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledVelocityTrajectory( ...
        ti, sampleData.uCenteredSubmesoscaleObserved(sampleIndices), sampleData.vCenteredSubmesoscaleObserved(sampleIndices), componentS, 0, 0, ...
        cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints);

    sampleStartIndex = sampleStartIndex + nSamples;
end

decomposition = struct( ...
    "fixedFrame", struct( ...
        "background", backgroundTrajectories, ...
        "mesoscale", mesoscaleTrajectories, ...
        "submesoscale", submesoscaleTrajectories), ...
    "centeredFrame", struct( ...
        "mesoscale", centeredMesoscaleTrajectories, ...
        "submesoscale", centeredSubmesoscaleTrajectories));
end

function [trajectory, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledVelocityTrajectory( ...
        t, uSamples, vSamples, componentS, xMean, yMean, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints)
[basisMatrix, knotPoints, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledSplineBasis( ...
    t, componentS, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints);
velocityCoefficients = basisMatrix \ [reshape(uSamples, [], 1), reshape(vSamples, [], 1)];
[xPathSpline, yPathSpline] = pathSplinesFromVelocityCoefficients(knotPoints, velocityCoefficients, componentS, xMean, yMean);
trajectory = TrajectorySpline(t=t, x=xPathSpline, y=yPathSpline);
end

function trajectory = reanchoredBackgroundTrajectory(backgroundTrajectory, ti)
xSpline = TensorSpline.fromKnotPoints(backgroundTrajectory.x.knotPoints, backgroundTrajectory.x.xi, ...
    S=backgroundTrajectory.x.S, xMean=-backgroundTrajectory.x(ti(1)));
ySpline = TensorSpline.fromKnotPoints(backgroundTrajectory.y.knotPoints, backgroundTrajectory.y.xi, ...
    S=backgroundTrajectory.y.S, xMean=-backgroundTrajectory.y(ti(1)));
trajectory = TrajectorySpline(t=ti, x=xSpline, y=ySpline);
end

function [basisMatrix, knotPoints, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints] = sampledSplineBasis( ...
        t, componentS, cachedTimes, cachedDegrees, cachedBasisMatrices, cachedKnotPoints)
t = reshape(t, [], 1);
for iCache = 1:numel(cachedTimes)
    if cachedDegrees(iCache) == componentS && isequal(cachedTimes{iCache}, t)
        basisMatrix = cachedBasisMatrices{iCache};
        knotPoints = cachedKnotPoints{iCache};
        return
    end
end

knotPoints = BSpline.knotPointsForDataPoints(t, S=componentS);
basisMatrix = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=componentS);
cachedTimes{end + 1, 1} = t;
cachedDegrees(end + 1, 1) = componentS;
cachedBasisMatrices{end + 1, 1} = basisMatrix;
cachedKnotPoints{end + 1, 1} = knotPoints;
end

function [xPathSpline, yPathSpline] = pathSplinesFromVelocityCoefficients(knotPoints, velocityCoefficients, splineDegree, xMean, yMean)
[integratedCoefficients, integratedKnotPoints, integratedS] = BSpline.integratedSplineState( ...
    velocityCoefficients, knotPoints=knotPoints, S=splineDegree);
xPathSpline = TensorSpline.fromKnotPoints(integratedKnotPoints, integratedCoefficients(:, 1), S=integratedS, xMean=xMean);
yPathSpline = TensorSpline.fromKnotPoints(integratedKnotPoints, integratedCoefficients(:, 2), S=integratedS, xMean=yMean);
end
