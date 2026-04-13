function diagnostics = referenceBandwidthForTwoDimensionalData(data, options)
% referenceBandwidthForTwoDimensionalData Solve the 2-D Botev bandwidth equations.

arguments
    data (:,2) double {mustBeReal, mustBeFinite, mustBeNonempty}
    options.minimum double {mustBeReal, mustBeFinite} = []
    options.maximum double {mustBeReal, mustBeFinite} = []
    options.bandwidthGridSize = []
end

sampleCount = size(data, 1);
if sampleCount <= size(data, 2)
    error("referenceBandwidthForTwoDimensionalData:TooFewSamples", ...
        "Two-dimensional density estimation requires more samples than dimensions.");
end

[minimum, maximum] = resolvedBoundsReference(data, options.minimum, options.maximum);
gridSize = resolvedBandwidthGridSizeReference(options.bandwidthGridSize);
scaling = maximum - minimum;
transformedData = (data - minimum)./scaling;
initialData = normalizedHistogram2DReference(transformedData, gridSize);
aSquared = discreteCosineTransform2DReference(initialData).^2;
iSquared = (0:gridSize - 1).^2;
halfWeights = [1 0.5 * ones(1, gridSize - 1)];
iPowers = cell(1, 6);
for iOrder = 0:5
    iPowers{iOrder + 1} = iSquared.^iOrder;
end

funcCacheKeys = cell(6, 6);
funcCacheValues = cell(6, 6);
psiCacheKeys = cell(6, 6);
psiCacheValues = cell(6, 6);
tStar = smallestRootReference(@objective, sampleCount);
p02 = func([0 2], tStar);
p20 = func([2 0], tStar);
p11 = func([1 1], tStar);
tY = (p02^(3/4)/(4*pi*sampleCount*p20^(3/4)*(p11 + sqrt(p20*p02))))^(1/3);
tX = (p20^(3/4)/(4*pi*sampleCount*p02^(3/4)*(p11 + sqrt(p20*p02))))^(1/3);
bandwidth = sqrt([tX tY]).*scaling;
objectiveResidual = objective(tStar);

diagnostics = struct( ...
    "bandwidth", bandwidth, ...
    "tStar", tStar, ...
    "p02", p02, ...
    "p20", p20, ...
    "p11", p11, ...
    "objectiveResidual", objectiveResidual, ...
    "minimum", minimum, ...
    "maximum", maximum, ...
    "bandwidthGridSize", gridSize);

    function value = objective(t)
        value = t - evolveResidual(t);
    end

    function value = evolveResidual(t)
        sumFunc = func([0 2], t) + func([2 0], t) + 2 * func([1 1], t);
        timeValue = (2*pi*sampleCount*sumFunc)^(-1/3);
        value = (t - timeValue)/timeValue;
    end

    function value = func(orderVector, t)
        iOrderX = orderVector(1) + 1;
        iOrderY = orderVector(2) + 1;
        cacheKey = sprintf("%.17g", t);
        [isCached, value] = lookupCachedValue(funcCacheKeys{iOrderX, iOrderY}, funcCacheValues{iOrderX, iOrderY}, cacheKey);
        if isCached
            return
        end

        if sum(orderVector) <= 4
            sumFunc = func([orderVector(1) + 1, orderVector(2)], t) + func([orderVector(1), orderVector(2) + 1], t);
            const = (1 + 1/2^(sum(orderVector) + 1))/3;
            timeValue = (-2*const*kernelConstantReference(orderVector(1))*kernelConstantReference(orderVector(2))/sampleCount/sumFunc)^(1/(2 + sum(orderVector)));
            value = psi(orderVector, timeValue);
        else
            value = psi(orderVector, t);
        end

        [funcCacheKeys{iOrderX, iOrderY}, funcCacheValues{iOrderX, iOrderY}] = appendCachedValue( ...
            funcCacheKeys{iOrderX, iOrderY}, funcCacheValues{iOrderX, iOrderY}, cacheKey, value);
    end

    function value = psi(orderVector, timeValue)
        iOrderX = orderVector(1) + 1;
        iOrderY = orderVector(2) + 1;
        cacheKey = sprintf("%.17g", timeValue);
        [isCached, value] = lookupCachedValue(psiCacheKeys{iOrderX, iOrderY}, psiCacheValues{iOrderX, iOrderY}, cacheKey);
        if isCached
            return
        end

        w = exp(-iSquared * pi^2 * timeValue) .* halfWeights;
        wx = w .* iPowers{orderVector(1) + 1};
        wy = w .* iPowers{orderVector(2) + 1};
        value = (-1)^sum(orderVector) * (wy * aSquared * wx.') * pi^(2 * sum(orderVector));
        [psiCacheKeys{iOrderX, iOrderY}, psiCacheValues{iOrderX, iOrderY}] = appendCachedValue( ...
            psiCacheKeys{iOrderX, iOrderY}, psiCacheValues{iOrderX, iOrderY}, cacheKey, value);
    end
end

function [minimum, maximum] = resolvedBoundsReference(data, minimumOption, maximumOption)
defaultMinimum = min(data, [], 1);
defaultMaximum = max(data, [], 1);
defaultRange = defaultMaximum - defaultMinimum;
defaultMinimum = defaultMinimum - defaultRange/2;
defaultMaximum = defaultMaximum + defaultRange/2;

if isempty(minimumOption)
    minimum = defaultMinimum;
else
    minimum = validatedBoundVectorReference(minimumOption, "minimum");
end

if isempty(maximumOption)
    maximum = defaultMaximum;
else
    maximum = validatedBoundVectorReference(maximumOption, "maximum");
end

if any(maximum <= minimum)
    error("referenceBandwidthForTwoDimensionalData:InvalidBounds", ...
        "maximum must be strictly greater than minimum in every dimension.");
end
end

function bound = validatedBoundVectorReference(bound, boundName)
bound = reshape(bound, 1, []);
if numel(bound) ~= 2
    error("referenceBandwidthForTwoDimensionalData:InvalidBounds", ...
        "%s must have exactly two entries.", boundName);
end
end

function gridSize = resolvedBandwidthGridSizeReference(bandwidthGridSize)
if isempty(bandwidthGridSize)
    gridSize = 2^8;
    return
end

if ~(isnumeric(bandwidthGridSize) && isreal(bandwidthGridSize) && isfinite(bandwidthGridSize) && isscalar(bandwidthGridSize) && bandwidthGridSize >= 1 && bandwidthGridSize == round(bandwidthGridSize))
    error("referenceBandwidthForTwoDimensionalData:InvalidBandwidthGridSize", ...
        "bandwidthGridSize must be a positive integer scalar.");
end

gridSize = 2^ceil(log2(double(bandwidthGridSize)));
end

function histogramValues = normalizedHistogram2DReference(data, gridSize)
nRows = size(data, 1);
bins = zeros(nRows, size(data, 2));
for iDimension = 1:size(data, 2)
    [~, bins(:, iDimension)] = histc(data(:, iDimension), 0:1/gridSize:1, 1);
    bins(:, iDimension) = min(bins(:, iDimension), gridSize);
end

histogramValues = accumarray(bins(all(bins > 0, 2), :), 1/nRows, [gridSize gridSize]);
end

function transformed = discreteCosineTransform2DReference(data)
[nRows, nCols] = size(data);
if nRows ~= nCols
    error("referenceBandwidthForTwoDimensionalData:InvalidGridSize", ...
        "The two-dimensional Botev solver requires a square DCT grid.");
end

weight = [1; 2 * exp(-1i * (1:(nRows - 1))' * pi/(2*nRows))];
weight = weight(:, ones(1, nCols));
transformed = discreteCosineTransformColumns(discreteCosineTransformColumns(data).').';

    function output = discreteCosineTransformColumns(input)
        output = real(weight .* fft([input(1:2:end, :); input(end:-2:2, :)]));
    end
end

function value = kernelConstantReference(order)
value = (-1)^order * prod(1:2:(2*order - 1))/sqrt(2*pi);
end

function value = smallestRootReference(objective, sampleCount)
sampleCount = 50 * (sampleCount <= 50) + 1050 * (sampleCount >= 1050) + sampleCount * ((sampleCount < 1050) & (sampleCount > 50));
searchLimit = 1e-12 + 0.01 * (sampleCount - 50)/1000;
lowerBound = 0;
objectiveLower = objective(lowerBound);

while true
    objectiveUpper = objective(searchLimit);
    if hasRootBracketReference(objectiveLower, objectiveUpper)
        value = fzero(objective, [lowerBound searchLimit]);
        return
    end

    if searchLimit == 0.1
        value = fminbnd(@(x) abs(objective(x)), 0, 0.1);
        return
    end

    searchLimit = min(searchLimit * 2, 0.1);
end
end

function tf = hasRootBracketReference(lowerValue, upperValue)
tf = false;
if ~isfinite(lowerValue) || ~isfinite(upperValue)
    return
end
if lowerValue == 0 || upperValue == 0
    tf = true;
    return
end
tf = sign(lowerValue) ~= sign(upperValue);
end

function [isCached, value] = lookupCachedValue(keys, values, cacheKey)
if isempty(keys)
    isCached = false;
    value = [];
    return
end

iCached = find(strcmp(keys, cacheKey), 1);
isCached = ~isempty(iCached);
if isCached
    value = values(iCached);
else
    value = [];
end
end

function [keys, values] = appendCachedValue(keys, values, cacheKey, value)
keys{end + 1, 1} = cacheKey;
values(end + 1, 1) = value;
end
