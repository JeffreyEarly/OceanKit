function theta = alignPrincipalStrainAngleToReference( ...
        sigma_n, sigma_s, referenceSigma_n, referenceSigma_s, options)
arguments
    sigma_n {mustBeNumeric, mustBeReal}
    sigma_s {mustBeNumeric, mustBeReal}
    referenceSigma_n {mustBeNumeric, mustBeReal}
    referenceSigma_s {mustBeNumeric, mustBeReal}
    options.units {mustBeTextScalar, mustBeMember(options.units, ["degrees", "radians"])} = "degrees"
end

if ~isequal(size(sigma_n), size(sigma_s))
    error("GriddedStreamfunction:StrainAngleSizeMismatch", ...
        "sigma_n and sigma_s must have the same size.");
end

principalTheta = principalStrainAngle(sigma_n, sigma_s);
referenceTheta = principalStrainAngle(referenceSigma_n, referenceSigma_s);

if isvector(principalTheta)
    thetaMatrix = reshape(principalTheta, [], 1);
    referenceMatrix = expandReferenceTheta(referenceTheta, numel(thetaMatrix), 1);
    thetaMatrix = alignThetaMatrix(thetaMatrix, referenceMatrix);
    theta = reshape(thetaMatrix, size(principalTheta));
else
    nTime = size(principalTheta, 1);
    thetaMatrix = reshape(principalTheta, nTime, []);
    referenceMatrix = expandReferenceTheta(referenceTheta, nTime, size(thetaMatrix, 2));
    thetaMatrix = alignThetaMatrix(thetaMatrix, referenceMatrix);
    theta = reshape(thetaMatrix, size(principalTheta));
end

if string(options.units) == "degrees"
    theta = theta * 180/pi;
end
end

function theta = principalStrainAngle(sigma_n, sigma_s)
theta = mod(atan2(sigma_s, sigma_n)/2 + pi/4, pi/2) - pi/4;
theta(~(isfinite(sigma_n) & isfinite(sigma_s))) = NaN;
end

function referenceTheta = expandReferenceTheta(referenceTheta, nTime, nSeries)
if isempty(referenceTheta)
    error("GriddedStreamfunction:EmptyReferenceTheta", ...
        "referenceSigma_n and referenceSigma_s must define at least one reference angle.");
end

referenceTheta = reshape(referenceTheta, size(referenceTheta, 1), []);
if isequal(size(referenceTheta), [nTime, nSeries])
    return
end

if isequal(size(referenceTheta), [nTime, 1])
    referenceTheta = repmat(referenceTheta, 1, nSeries);
    return
end

if isequal(size(referenceTheta), [1, nSeries])
    referenceTheta = repmat(referenceTheta, nTime, 1);
    return
end

if isscalar(referenceTheta)
    referenceTheta = repmat(referenceTheta, nTime, nSeries);
    return
end

error("GriddedStreamfunction:ReferenceThetaSizeMismatch", ...
    "Reference strain angle must be scalar, [nTime 1], [1 nSeries], or [nTime nSeries].");
end

function theta = alignThetaMatrix(theta, referenceTheta)
halfTurn = pi/2;
for iTime = 1:size(theta, 1)
    isValid = isfinite(theta(iTime, :)) & isfinite(referenceTheta(iTime, :));
    branchShift = zeros(1, size(theta, 2));
    branchShift(isValid) = round((referenceTheta(iTime, isValid) - theta(iTime, isValid))/halfTurn);
    theta(iTime, isValid) = theta(iTime, isValid) + branchShift(isValid) * halfTurn;
    theta(iTime, ~isValid) = NaN;
end
end
