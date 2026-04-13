function diagnostics = diagnosticsFromFit(fit, totalDisplacementWeightsByTrajectory)
if nargin < 2
    totalDisplacementWeightsByTrajectory = {};
end

diagnostics = struct( ...
    "kappa", GriddedStreamfunctionBootstrap.kappaEstimateFromFit(fit, totalDisplacementWeightsByTrajectory), ...
    "coherence", NaN, ...
    "coherenceSpectrum", emptyCoherenceSpectrum());

if ~coherenceBackendAvailable()
    return
end

tCell = arrayfun(@(trajectory) reshape(trajectory.t, [], 1), reshape(fit.observedTrajectories, [], 1), 'UniformOutput', false);
tCoherence = coherenceEvaluationTimes(tCell);
if numel(tCoherence) < 2
    return
end

centeredFrame = fit.decomposition.centeredFrame;
nTrajectories = numel(centeredFrame.mesoscale);
cvMesoscale = complex(zeros(numel(tCoherence), nTrajectories));
cvSubmesoscale = complex(zeros(numel(tCoherence), nTrajectories));

for iTrajectory = 1:nTrajectories
    mesoscale = centeredFrame.mesoscale(iTrajectory);
    submesoscale = centeredFrame.submesoscale(iTrajectory);
    cvMesoscale(:, iTrajectory) = reshape(mesoscale.u(tCoherence), [], 1) + 1i * reshape(mesoscale.v(tCoherence), [], 1);
    cvSubmesoscale(:, iTrajectory) = reshape(submesoscale.u(tCoherence), [], 1) + 1i * reshape(submesoscale.v(tCoherence), [], 1);
end

dtCoherence = tCoherence(2) - tCoherence(1);
[psi, ~] = sleptap(size(cvMesoscale, 1));
[frequency, sxx, syy, sxy] = mspec(dtCoherence, cvMesoscale, cvSubmesoscale, psi, 'cyclic');
gamma = abs(sxy).^2 ./ (sxx .* syy);
meanCoherence = meanOverFinite(gamma);
frequency = reshape(frequency, [], 1);
meanCoherence = reshape(meanCoherence, [], 1);
lowerFrequencyHalf = meanCoherence(1:ceil(numel(meanCoherence) / 2));
finiteCoherence = isfinite(lowerFrequencyHalf);
if ~any(finiteCoherence)
    return
end

diagnostics.coherence = mean(lowerFrequencyHalf(finiteCoherence));
diagnostics.coherenceSpectrum = struct( ...
    "frequency", frequency, ...
    "coherence", meanCoherence);
end

function spectrum = emptyCoherenceSpectrum()
spectrum = struct("frequency", zeros(0, 1), "coherence", zeros(0, 1));
end

function tf = coherenceBackendAvailable()
persistent didWarnMissing

tf = exist('sleptap', 'file') ~= 0 && exist('mspec', 'file') ~= 0;
if ~tf && isempty(didWarnMissing)
    warning('GriddedStreamfunctionBootstrap:MissingCoherenceBackend', ...
        'jLab functions sleptap and mspec were not found. Coherence diagnostics will be left unavailable.');
    didWarnMissing = true;
end
end

function tCoherence = coherenceEvaluationTimes(tCell)
nTrajectories = numel(tCell);
tStart = -Inf;
tEnd = Inf;
dtValues = zeros(nTrajectories, 1);

for iTrajectory = 1:nTrajectories
    ti = reshape(tCell{iTrajectory}, [], 1);
    if numel(ti) < 2
        tCoherence = zeros(0, 1);
        return
    end

    dtValues(iTrajectory) = median(diff(ti));
    if ~isfinite(dtValues(iTrajectory)) || dtValues(iTrajectory) <= 0
        tCoherence = zeros(0, 1);
        return
    end

    tStart = max(tStart, ti(1));
    tEnd = min(tEnd, ti(end));
end

if ~(isfinite(tStart) && isfinite(tEnd)) || tEnd <= tStart
    tCoherence = zeros(0, 1);
    return
end

dtCoherence = max(dtValues);
nStep = floor((tEnd - tStart) / dtCoherence);
if nStep < 1
    tCoherence = zeros(0, 1);
    return
end

tCoherence = tStart + (0:nStep).' * dtCoherence;
end

function meanValues = meanOverFinite(values)
meanValues = zeros(size(values, 1), 1);
for iRow = 1:size(values, 1)
    finiteValues = values(iRow, isfinite(values(iRow, :)));
    if isempty(finiteValues)
        meanValues(iRow) = NaN;
    else
        meanValues(iRow) = mean(finiteValues);
    end
end
end
