function sampleData = decompositionSampleData(self, trajectories)
arguments (Input)
    self (1,1) GriddedStreamfunction
    trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
end
arguments (Output)
    sampleData struct
end

trajectories = reshape(trajectories, [], 1);
if shouldReuseObservedTrajectorySamples(self, trajectories)
    sampleData = self.observedTrajectorySampleData;
else
    sampleData = GriddedStreamfunction.sampleTrajectoryData(trajectories);
end

allT = sampleData.allT;
allX = sampleData.allX;
allY = sampleData.allY;
allXDot = sampleData.allXDot;
allYDot = sampleData.allYDot;

validateFastTimeDomain(self, allT);

mxDotAll = reshape(self.centerOfMassTrajectory.u(allT), [], 1);
myDotAll = reshape(self.centerOfMassTrajectory.v(allT), [], 1);
[~, qAll, rAll] = centeredCoordinates(self, allT, allX, allY);
validatePsiDomain(self, qAll, rAll, allT);

uMesoscaleObserved = -self.streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[0 1 0]);
vMesoscaleObserved = self.streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[1 0 0]);
uMesoscaleObserved = reshape(uMesoscaleObserved, [], 1);
vMesoscaleObserved = reshape(vMesoscaleObserved, [], 1);

fastState = GriddedStreamfunction.fastBasisState(allT, self.fastKnotPoints, self.fastS, false);
% Reuse the fast-basis factorization in the same least-squares projector B * (B \ A).
mesoscaleComObserved = fastState.B * (fastState.solver \ [uMesoscaleObserved, vMesoscaleObserved]);
uMesoscaleComObserved = mesoscaleComObserved(:, 1);
vMesoscaleComObserved = mesoscaleComObserved(:, 2);
uMesoscaleRelativeObserved = uMesoscaleObserved - uMesoscaleComObserved;
vMesoscaleRelativeObserved = vMesoscaleObserved - vMesoscaleComObserved;

uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
vBackgroundObserved = myDotAll - vMesoscaleComObserved;
uCenteredObserved = allXDot - mxDotAll;
vCenteredObserved = allYDot - myDotAll;
uCenteredSubmesoscaleObserved = uCenteredObserved - uMesoscaleRelativeObserved;
vCenteredSubmesoscaleObserved = vCenteredObserved - vMesoscaleRelativeObserved;
uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

sampleData.mxDotAll = mxDotAll;
sampleData.myDotAll = myDotAll;
sampleData.qAll = qAll;
sampleData.rAll = rAll;
sampleData.B = fastState.B;
sampleData.uMesoscaleObserved = uMesoscaleObserved;
sampleData.vMesoscaleObserved = vMesoscaleObserved;
sampleData.uMesoscaleComObserved = uMesoscaleComObserved;
sampleData.vMesoscaleComObserved = vMesoscaleComObserved;
sampleData.uMesoscaleRelativeObserved = uMesoscaleRelativeObserved;
sampleData.vMesoscaleRelativeObserved = vMesoscaleRelativeObserved;
sampleData.uBackgroundObserved = uBackgroundObserved;
sampleData.vBackgroundObserved = vBackgroundObserved;
sampleData.uCenteredObserved = uCenteredObserved;
sampleData.vCenteredObserved = vCenteredObserved;
sampleData.uCenteredSubmesoscaleObserved = uCenteredSubmesoscaleObserved;
sampleData.vCenteredSubmesoscaleObserved = vCenteredSubmesoscaleObserved;
sampleData.uSubmesoscaleObserved = uSubmesoscaleObserved;
sampleData.vSubmesoscaleObserved = vSubmesoscaleObserved;
end

function tf = shouldReuseObservedTrajectorySamples(self, trajectories)
if isempty(self.observedTrajectorySampleData) || numel(trajectories) ~= numel(self.observedTrajectories)
    tf = false;
    return
end

tf = true;
for iTrajectory = 1:numel(trajectories)
    if trajectories(iTrajectory) ~= self.observedTrajectories(iTrajectory)
        tf = false;
        return
    end
end
end

function validateFastTimeDomain(self, allT)
if any(allT < self.fastKnotPoints(1) | allT > self.fastKnotPoints(end))
    error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
        "Observation times must lie inside the fitted fast spline domain.");
end
end

function validatePsiDomain(self, qAll, rAll, allT)
qDomain = [self.psiKnotPoints{1}(1), self.psiKnotPoints{1}(end)];
rDomain = [self.psiKnotPoints{2}(1), self.psiKnotPoints{2}(end)];
tDomain = [self.psiKnotPoints{3}(1), self.psiKnotPoints{3}(end)];

if any(qAll < qDomain(1) | qAll > qDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered q observations must lie inside the fitted psi spline domain.");
end
if any(rAll < rDomain(1) | rAll > rDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered r observations must lie inside the fitted psi spline domain.");
end
if any(allT < tDomain(1) | allT > tDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Observation times must lie inside the fitted psi time domain.");
end
end
