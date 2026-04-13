function scores = computeComparisonConsensusScores(summary, scoreIndices, mesoscaleConstraint)
arguments
    summary struct
    scoreIndices {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    mesoscaleConstraint {mustBeTextScalar, mustBeMember(mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
end

requiredFields = ["uCenter", "vCenter", "sigma_n", "sigma_s"];
for iField = 1:numel(requiredFields)
    if ~isfield(summary, requiredFields(iField))
        error("GriddedStreamfunction:MissingComparisonSummaryField", ...
            "summary must contain the field '%s'.", requiredFields(iField));
    end
end

if mesoscaleConstraint ~= "zeroVorticity" && ~isfield(summary, "zeta")
    error("GriddedStreamfunction:MissingComparisonSummaryField", ...
        "summary must contain the field 'zeta' when mesoscaleConstraint is '%s'.", mesoscaleConstraint);
end

scoreIndices = reshape(scoreIndices, 1, []);
nTime = size(summary.uCenter, 1);
if any(scoreIndices < 1 | scoreIndices > nTime | scoreIndices ~= round(scoreIndices))
    error("GriddedStreamfunction:InvalidScoreIndices", ...
        "scoreIndices must be integer indices between 1 and %d.", nTime);
end

nBootstraps = size(summary.uCenter, 2);
scores = struct( ...
    "uv", zeros(1, nBootstraps), ...
    "strain", zeros(1, nBootstraps), ...
    "zeta", zeros(1, nBootstraps), ...
    "joint", zeros(1, nBootstraps));

for iScore = scoreIndices
    scores.uv = scores.uv + score2DBlock(summary.uCenter(iScore, :), summary.vCenter(iScore, :));

    if mesoscaleConstraint ~= "zeroStrain"
        scores.strain = scores.strain + score2DBlock(summary.sigma_n(iScore, :), summary.sigma_s(iScore, :));
    end

    if mesoscaleConstraint ~= "zeroVorticity"
        scores.zeta = scores.zeta + score1DBlock(summary.zeta(iScore, :));
    end
end

scores.joint = scores.uv + scores.strain + scores.zeta;
end

function score = score2DBlock(xValues, yValues)
xValues = reshape(xValues, [], 1);
yValues = reshape(yValues, [], 1);
score = zeros(1, numel(xValues));

[isConstantX, isConstantY] = constantFlags(xValues, yValues);
if isConstantX && isConstantY
    return
end

if isConstantX
    score = score1DBlock(yValues);
    return
end

if isConstantY
    score = score1DBlock(xValues);
    return
end

model = KernelDensityEstimate.fromData([xValues, yValues]);
pointDensity = model.densityAt([xValues, yValues]);
score = log10(safeDensity(pointDensity)).';
end

function score = score1DBlock(values)
values = reshape(values, [], 1);
score = zeros(1, numel(values));

if isConstant(values)
    return
end

model = KernelDensityEstimate.fromData(values);
pointDensity = model.densityAt(values);
score = log10(safeDensity(pointDensity)).';
end

function density = safeDensity(values)
density = reshape(values, [], 1);
density(~isfinite(density) | density <= 0) = realmin("double");
end

function [isConstantX, isConstantY] = constantFlags(xValues, yValues)
isConstantX = isConstant(xValues);
isConstantY = isConstant(yValues);
end

function tf = isConstant(values)
values = reshape(values, [], 1);
scale = max([1; abs(values)]);
tf = (max(values) - min(values)) <= 1e-12 * scale;
end
