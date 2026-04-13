classdef KernelDensityEstimate
    % Fit and evaluate Gaussian kernel density estimates in one or two dimensions.
    %
    % `KernelDensityEstimate` stores sampled data together with a Gaussian
    % product-kernel bandwidth selected from the diffusion-based Botev
    % estimator. For fitted models created with `fromData(...)`, the
    % bandwidth is chosen from the same diffusion fixed-point equations as
    % the legacy `kde` and `kde2d` utilities, while density and CDF values
    % are evaluated exactly from the Gaussian kernel sum rather than from a
    % reconstructed DCT grid.
    %
    % In one dimension the fitted density is
    %
    % $$ \hat{f}(x) = \frac{1}{N h}\sum_{i=1}^{N}\phi\!\left(\frac{x - x_i}{h}\right), $$
    %
    % and in two dimensions the fitted density is
    %
    % $$ \hat{f}(x,y) = \frac{1}{N h_x h_y}\sum_{i=1}^{N}\phi\!\left(\frac{x - x_i}{h_x}\right)\phi\!\left(\frac{y - y_i}{h_y}\right), $$
    %
    % where $$\phi$$ is the standard normal density.
    %
    % ```matlab
    % data = [sigma_n(:), sigma_s(:)];
    % model = KernelDensityEstimate.fromData(data);
    % [density, gridVectors] = model.densityOnGrid(gridSize=[192 320]);
    % contourf(gridVectors{1}, gridVectors{2}, density.')
    % stats = KernelDensityEstimate.planarStatisticsFromData(data);
    % polar = KernelDensityEstimate.polarSummaryFromPlanarStatistics(stats);
    % strain = KernelDensityEstimate.strainSummaryFromPlanarStatistics(stats);
    % ```
    %
    % - Topic: Create a density estimate
    % - Topic: Inspect density-estimate properties
    % - Topic: Evaluate the density estimate
    % - Topic: Evaluate the density estimate — Planar summary
    % - Topic: Evaluate the density estimate — Planar rendering
    % - Topic: Evaluate the density estimate — Polar reduction
    % - Topic: Evaluate the density estimate — Strain reduction
    % - Declaration: classdef KernelDensityEstimate

    properties (SetAccess = private)
        % Sample locations used to define the fitted kernel density estimate.
        %
        % `data` is stored as an `N-by-1` column vector for one-dimensional
        % fits and as an `N-by-2` matrix for two-dimensional fits.
        %
        % - Topic: Inspect density-estimate properties
        data (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty} = 0

        % Gaussian kernel bandwidth.
        %
        % `bandwidth` is scalar in one dimension and a `1-by-2` row vector
        % `[hx hy]` in two dimensions.
        %
        % - Topic: Inspect density-estimate properties
        bandwidth (1,:) double {mustBeReal, mustBeFinite, mustBePositive} = 1

        % Lower bound of the default evaluation box.
        %
        % `minimum` stores the lower endpoint used by `densityOnGrid(...)`
        % and `cdfOnGrid(...)`.
        %
        % - Topic: Inspect density-estimate properties
        minimum (1,:) double {mustBeReal, mustBeFinite} = 0

        % Upper bound of the default evaluation box.
        %
        % `maximum` stores the upper endpoint used by `densityOnGrid(...)`
        % and `cdfOnGrid(...)`.
        %
        % - Topic: Inspect density-estimate properties
        maximum (1,:) double {mustBeReal, mustBeFinite} = 1
    end

    properties (Dependent, SetAccess = private)
        % Number of coordinates per sample.
        %
        % `numDimensions` is `1` for vector-valued sample data and `2` for
        % planar sample clouds.
        %
        % - Topic: Inspect density-estimate properties
        numDimensions
    end

    methods
        function self = KernelDensityEstimate(options)
            % Create a kernel density estimate from solved state.
            %
            % Use this cheap canonical constructor when the sample data,
            % bandwidth, and default evaluation box are already known.
            % Ordinary fitting workflows should prefer `fromData(...)`.
            %
            % - Topic: Create a density estimate
            % - Declaration: self = KernelDensityEstimate(data=...,bandwidth=...,minimum=...,maximum=...)
            % - Parameter data: sample locations stored as an `N-by-1` vector or `N-by-2` matrix
            % - Parameter bandwidth: scalar Gaussian bandwidth in one dimension or `[hx hy]` in two dimensions
            % - Parameter minimum: lower endpoint of the default evaluation box
            % - Parameter maximum: upper endpoint of the default evaluation box
            % - Returns self: canonical `KernelDensityEstimate` instance
            arguments
                options.data (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty}
                options.bandwidth (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
                options.minimum (1,:) double {mustBeReal, mustBeFinite}
                options.maximum (1,:) double {mustBeReal, mustBeFinite}
            end

            data = normalizedSampleData(options.data);
            numDimensions = size(data, 2);
            validateSupportedDimensionCount(numDimensions);
            minimum = validatedBoundVector(options.minimum, numDimensions, "minimum");
            maximum = validatedBoundVector(options.maximum, numDimensions, "maximum");
            validateStrictlyIncreasingBounds(minimum, maximum);
            validateBandwidthVector(options.bandwidth, numDimensions);

            self.data = data;
            self.bandwidth = options.bandwidth;
            self.minimum = minimum;
            self.maximum = maximum;
        end

        function value = get.numDimensions(self)
            value = size(self.data, 2);
        end

        function density = densityAt(self, queryPoints)
            % Evaluate the fitted density at arbitrary query points.
            %
            % In one dimension `queryPoints` may be any numeric array, and
            % the returned `density` matches its shape. In two dimensions
            % `queryPoints` must be an `Nq-by-2` array and the returned
            % `density` is an `Nq-by-1` column vector.
            %
            % - Topic: Evaluate the density estimate
            % - Declaration: density = densityAt(self,queryPoints)
            % - Parameter queryPoints: evaluation points as a numeric array in one dimension or an `Nq-by-2` array in two dimensions
            % - Returns density: exact Gaussian KDE values at the requested query points
            arguments
                self (1,1) KernelDensityEstimate
                queryPoints {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonempty}
            end

            [queryMatrix, queryShape] = normalizedQueryPoints(queryPoints, self.numDimensions);

            if self.numDimensions == 1
                density = gaussianDensityAt1D(self.data, self.bandwidth, queryMatrix);
                density = reshape(density, queryShape);
            else
                density = gaussianDensityAt2D(self.data, self.bandwidth, queryMatrix);
            end
        end

        function [density, gridVectors] = densityOnGrid(self, options)
            % Evaluate the fitted density on a regular query grid.
            %
            % `densityOnGrid(...)` uses the stored default bounds together
            % with an arbitrary positive evaluation grid size. The public
            % evaluation grid does not need to be a power of two.
            %
            % - Topic: Evaluate the density estimate
            % - Declaration: [density,gridVectors] = densityOnGrid(self,gridSize=...)
            % - Parameter gridSize: optional number of grid points; scalar in one dimension, or scalar or `[Nx Ny]` in two dimensions
            % - Returns density: exact Gaussian KDE values on the regular query grid
            % - Returns gridVectors: cell array of one or two column vectors defining the regular query grid
            arguments
                self (1,1) KernelDensityEstimate
                options.gridSize = []
            end

            gridVectors = evaluationGridVectors(self.minimum, self.maximum, self.numDimensions, options.gridSize);

            if self.numDimensions == 1
                density = self.densityAt(gridVectors{1});
                return
            end

            [X, Y] = ndgrid(gridVectors{1}, gridVectors{2});
            density = reshape(self.densityAt([X(:), Y(:)]), size(X));
        end

        function [cdfValues, gridVectors] = cdfOnGrid(self, options)
            % Evaluate the fitted one-dimensional CDF on a regular query grid.
            %
            % The returned CDF is the exact Gaussian-mixture cumulative
            % distribution associated with the stored sample locations and
            % bandwidth.
            %
            % - Topic: Evaluate the density estimate
            % - Declaration: [cdfValues,gridVectors] = cdfOnGrid(self,gridSize=...)
            % - Parameter gridSize: optional number of grid points in the regular one-dimensional query grid
            % - Returns cdfValues: exact Gaussian-mixture CDF values on the regular query grid
            % - Returns gridVectors: one-element cell array containing the query grid as a column vector
            arguments
                self (1,1) KernelDensityEstimate
                options.gridSize = []
            end

            if self.numDimensions ~= 1
                error("KernelDensityEstimate:InvalidDimension", ...
                    "cdfOnGrid is only defined for one-dimensional kernel density estimates.");
            end

            gridVectors = evaluationGridVectors(self.minimum, self.maximum, self.numDimensions, options.gridSize);
            cdfValues = gaussianCdfAt1D(self.data, self.bandwidth, gridVectors{1});
        end
    end

    methods (Static)
        function self = fromData(data, options)
            % Fit a kernel density estimate from sampled data.
            %
            % `fromData(...)` chooses a Gaussian product-kernel bandwidth
            % using the diffusion estimator of Botev, Grotowski, and
            % Kroese, then stores the fitted data cloud together with the
            % resulting bandwidth and default evaluation box.
            %
            % - Topic: Create a density estimate
            % - nav_order: 1
            % - Declaration: self = KernelDensityEstimate.fromData(data,minimum=...,maximum=...,bandwidthGridSize=...)
            % - Parameter data: sample locations stored as a vector for one dimension or an `N-by-2` matrix for two dimensions
            % - Parameter minimum: optional lower endpoint of the default evaluation box; the default is `min(data) - range(data)/2` componentwise
            % - Parameter maximum: optional upper endpoint of the default evaluation box; the default is `max(data) + range(data)/2` componentwise
            % - Parameter bandwidthGridSize: optional DCT grid size used internally by the Botev solver; values are rounded up internally to the next power of two
            % - Returns self: fitted `KernelDensityEstimate` instance
            arguments
                data (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty}
                options.minimum double {mustBeReal, mustBeFinite} = []
                options.maximum double {mustBeReal, mustBeFinite} = []
                options.bandwidthGridSize = []
            end

            data = normalizedSampleData(data);
            numDimensions = size(data, 2);
            validateSupportedDimensionCount(numDimensions);
            [minimum, maximum] = resolvedBounds(data, options.minimum, options.maximum);
            gridSize = resolvedBandwidthGridSize(numDimensions, options.bandwidthGridSize);

            if numDimensions == 1
                bandwidth = bandwidthForOneDimensionalData(data, minimum, maximum, gridSize);
            else
                bandwidth = bandwidthForTwoDimensionalData(data, minimum, maximum, gridSize);
            end

            self = KernelDensityEstimate(data=data, bandwidth=bandwidth, minimum=minimum, maximum=maximum);
        end

        statistics = planarStatisticsFromData(data, options)
        polarSummary = polarSummaryFromPlanarStatistics(statistics, options)
        strainSummary = strainSummaryFromPlanarStatistics(statistics, options)
        plotPlanarStatistics(ax, statistics, options)
    end
end

function data = normalizedSampleData(data)
if isvector(data)
    data = reshape(data, [], 1);
end
end

function validateSupportedDimensionCount(numDimensions)
if ~ismember(numDimensions, [1 2])
    error("KernelDensityEstimate:UnsupportedDimension", ...
        "data must be a vector of one-dimensional samples or an N-by-2 array of two-dimensional samples.");
end
end

function bound = validatedBoundVector(bound, numDimensions, boundName)
bound = reshape(bound, 1, []);
if numel(bound) ~= numDimensions
    error("KernelDensityEstimate:InvalidBounds", ...
        "%s must have one entry per sample dimension.", boundName);
end
end

function validateStrictlyIncreasingBounds(minimum, maximum)
if any(maximum <= minimum)
    error("KernelDensityEstimate:InvalidBounds", ...
        "maximum must be strictly greater than minimum in every dimension.");
end
end

function validateBandwidthVector(bandwidth, numDimensions)
if numel(bandwidth) ~= numDimensions
    error("KernelDensityEstimate:InvalidBandwidth", ...
        "bandwidth must have one entry per sample dimension.");
end
end

function [minimum, maximum] = resolvedBounds(data, minimumOption, maximumOption)
numDimensions = size(data, 2);
defaultMinimum = min(data, [], 1);
defaultMaximum = max(data, [], 1);
defaultRange = defaultMaximum - defaultMinimum;
defaultMinimum = defaultMinimum - defaultRange/2;
defaultMaximum = defaultMaximum + defaultRange/2;

if isempty(minimumOption)
    minimum = defaultMinimum;
else
    minimum = validatedBoundVector(minimumOption, numDimensions, "minimum");
end

if isempty(maximumOption)
    maximum = defaultMaximum;
else
    maximum = validatedBoundVector(maximumOption, numDimensions, "maximum");
end

validateStrictlyIncreasingBounds(minimum, maximum)
end

function gridSize = resolvedBandwidthGridSize(numDimensions, bandwidthGridSize)
if isempty(bandwidthGridSize)
    if numDimensions == 1
        gridSize = 2^14;
    else
        gridSize = 2^8;
    end
    return
end

if ~(isnumeric(bandwidthGridSize) && isreal(bandwidthGridSize) && isfinite(bandwidthGridSize) && isscalar(bandwidthGridSize) && bandwidthGridSize >= 1 && bandwidthGridSize == round(bandwidthGridSize))
    error("KernelDensityEstimate:InvalidBandwidthGridSize", ...
        "bandwidthGridSize must be a positive integer scalar.");
end

gridSize = roundUpPowerOfTwo(bandwidthGridSize);
end

function gridVectors = evaluationGridVectors(minimum, maximum, numDimensions, gridSizeOption)
gridSize = resolvedEvaluationGridSize(numDimensions, gridSizeOption);
gridVectors = cell(1, numDimensions);
for iDimension = 1:numDimensions
    gridVectors{iDimension} = linspace(minimum(iDimension), maximum(iDimension), gridSize(iDimension)).';
end
end

function gridSize = resolvedEvaluationGridSize(numDimensions, gridSizeOption)
if isempty(gridSizeOption)
    if numDimensions == 1
        gridSize = 16384;
    else
        gridSize = [256 256];
    end
elseif numDimensions == 1
    gridSize = validateGridSizeScalar(gridSizeOption);
else
    gridSize = validateGridSizeVector(gridSizeOption);
end
end

function gridSize = validateGridSizeScalar(gridSize)
if ~(isnumeric(gridSize) && isreal(gridSize) && isfinite(gridSize) && isscalar(gridSize) && gridSize >= 1 && gridSize == round(gridSize))
    error("KernelDensityEstimate:InvalidGridSize", ...
        "gridSize must be a positive integer.");
end

gridSize = double(gridSize);
end

function gridSize = validateGridSizeVector(gridSize)
if ~(isnumeric(gridSize) && isreal(gridSize) && all(isfinite(gridSize), "all"))
    error("KernelDensityEstimate:InvalidGridSize", ...
        "gridSize must be a positive integer scalar or a 1-by-2 vector of positive integers.");
end

gridSize = reshape(double(gridSize), 1, []);
if isscalar(gridSize)
    gridSize = [gridSize gridSize];
elseif numel(gridSize) ~= 2
    error("KernelDensityEstimate:InvalidGridSize", ...
        "gridSize must be a positive integer scalar or a 1-by-2 vector of positive integers.");
end

if any(gridSize < 1 | gridSize ~= round(gridSize))
    error("KernelDensityEstimate:InvalidGridSize", ...
        "gridSize must be a positive integer scalar or a 1-by-2 vector of positive integers.");
end
end

function [queryMatrix, queryShape] = normalizedQueryPoints(queryPoints, numDimensions)
if numDimensions == 1
    queryShape = size(queryPoints);
    queryMatrix = reshape(queryPoints, [], 1);
    return
end

queryShape = [numel(queryPoints) 1];
if isvector(queryPoints) && numel(queryPoints) == 2
    queryMatrix = reshape(queryPoints, 1, 2);
    return
end

if size(queryPoints, 2) ~= 2
    error("KernelDensityEstimate:InvalidQueryPoints", ...
        "Two-dimensional queryPoints must be an N-by-2 array.");
end

queryMatrix = queryPoints;
end

function density = gaussianDensityAt1D(data, bandwidth, queryPoints)
nSamples = size(data, 1);
queryPoints = reshape(queryPoints, [], 1);
density = zeros(size(queryPoints));
normalization = bandwidth * sqrt(2*pi);
blockSize = queryBlockSize(nSamples, 1);

for iStart = 1:blockSize:numel(queryPoints)
    iEnd = min(iStart + blockSize - 1, numel(queryPoints));
    scaledOffset = (queryPoints(iStart:iEnd) - data.')/bandwidth;
    density(iStart:iEnd) = mean(exp(-0.5 * scaledOffset.^2), 2)/normalization;
end
end

function density = gaussianDensityAt2D(data, bandwidth, queryPoints)
nSamples = size(data, 1);
density = zeros(size(queryPoints, 1), 1);
normalization = 2*pi * bandwidth(1) * bandwidth(2);
blockSize = queryBlockSize(nSamples, 2);

for iStart = 1:blockSize:size(queryPoints, 1)
    iEnd = min(iStart + blockSize - 1, size(queryPoints, 1));
    scaledX = (queryPoints(iStart:iEnd, 1) - data(:, 1).')/bandwidth(1);
    scaledY = (queryPoints(iStart:iEnd, 2) - data(:, 2).')/bandwidth(2);
    density(iStart:iEnd) = mean(exp(-0.5 * (scaledX.^2 + scaledY.^2)), 2)/normalization;
end
end

function cdfValues = gaussianCdfAt1D(data, bandwidth, queryPoints)
nSamples = size(data, 1);
queryPoints = reshape(queryPoints, [], 1);
cdfValues = zeros(size(queryPoints));
blockSize = queryBlockSize(nSamples, 1);
normalization = bandwidth * sqrt(2);

for iStart = 1:blockSize:numel(queryPoints)
    iEnd = min(iStart + blockSize - 1, numel(queryPoints));
    scaledOffset = (queryPoints(iStart:iEnd) - data.')/normalization;
    cdfValues(iStart:iEnd) = mean(0.5 * erfc(-scaledOffset), 2);
end
end

function blockSize = queryBlockSize(nSamples, numDimensions)
targetElements = 4e6;
blockSize = max(1, floor(targetElements/(nSamples * numDimensions)));
end

function bandwidth = bandwidthForOneDimensionalData(data, minimum, maximum, gridSize)
data = reshape(data, [], 1);
rangeWidth = maximum - minimum;
xMesh = linspace(minimum, maximum, gridSize).';
sampleCount = length(unique(data));
initialData = histc(data, xMesh)/sampleCount;
initialData = initialData/sum(initialData);
a = discreteCosineTransform1D(initialData);
iSquared = (1:(gridSize - 1))'.^2;
aSquared = (a(2:end)/2).^2;
maxOrder = 7;
momentOrders = 2:maxOrder;
expScale = iSquared * pi^2;
orderMomentWeights = zeros(maxOrder - 1, numel(iSquared));
for iOrder = 1:length(momentOrders)
    orderMomentWeights(iOrder, :) = 2*pi^(2*momentOrders(iOrder)) * ((iSquared.^momentOrders(iOrder)) .* aSquared).';
end
recursiveNumerators = zeros(1, maxOrder - 2);
recursiveExponents = zeros(1, maxOrder - 2);
for s = 2:(maxOrder - 1)
    k0 = prod(1:2:(2*s - 1))/sqrt(2*pi);
    plugInConstant = (1 + (1/2)^(s + 1/2))/3;
    recursiveNumerators(s - 1) = 2*plugInConstant*k0/sampleCount;
    recursiveExponents(s - 1) = 2/(3 + 2*s);
end
finalResidualScale = 2*sampleCount*sqrt(pi);
cachedTimes = zeros(0, 1);
cachedResiduals = zeros(0, 1);
tStar = smallestRoot(@fixedPointResidual, sampleCount);
bandwidth = sqrt(tStar) * rangeWidth;

    function residual = fixedPointResidual(t)
        iCached = find(cachedTimes == t, 1);
        if ~isempty(iCached)
            residual = cachedResiduals(iCached);
            return
        end

        momentValue = derivativeMoment(maxOrder, t);
        for s = (maxOrder - 1):-1:2
            timeValue = (recursiveNumerators(s - 1)/momentValue)^recursiveExponents(s - 1);
            momentValue = derivativeMoment(s, timeValue);
        end
        residual = t - (finalResidualScale*momentValue)^(-2/5);
        cachedTimes(end + 1, 1) = t; %#ok<AGROW>
        cachedResiduals(end + 1, 1) = residual; %#ok<AGROW>
    end

    function momentValue = derivativeMoment(order, timeValue)
        momentValue = orderMomentWeights(order - 1, :) * exp(-expScale * timeValue);
    end
end

function bandwidth = bandwidthForTwoDimensionalData(data, minimum, maximum, gridSize)
sampleCount = size(data, 1);
if sampleCount <= size(data, 2)
    error("KernelDensityEstimate:TooFewSamples", ...
        "Two-dimensional density estimation requires more samples than dimensions.");
end

scaling = maximum - minimum;
transformedData = (data - minimum)./scaling;
initialData = normalizedHistogram2D(transformedData, gridSize);
a = discreteCosineTransform2D(initialData);
aSquared = a.^2;
iSquared = (0:gridSize - 1).^2;
halfWeights = [1 0.5 * ones(1, gridSize - 1)];
iPowers = cell(1, 6);
for iOrder = 0:5
    iPowers{iOrder + 1} = iSquared.^iOrder;
end
piSquared = pi^2;
weightedSpectrum = (halfWeights.' * halfWeights) .* aSquared;
orderPowerMatrix = vertcat(iPowers{:});
orderConstants = zeros(1, 5);
for sumOrder = 0:4
    orderConstants(sumOrder + 1) = (1 + 1/2^(sumOrder + 1))/3;
end
kernelConstants = zeros(1, 6);
for iOrder = 0:5
    kernelConstants(iOrder + 1) = kernelConstant(iOrder);
end
psiScale = zeros(6, 6);
for iOrderY = 0:5
    for iOrderX = 0:5
        psiScale(iOrderY + 1, iOrderX + 1) = (-1)^(iOrderX + iOrderY) * pi^(2 * (iOrderX + iOrderY));
    end
end
funcTimeNumerator = zeros(6, 6);
funcTimeExponent = zeros(6, 6);
for iOrderY = 0:5
    for iOrderX = 0:5
        sumOrder = iOrderX + iOrderY;
        if sumOrder <= 4
            funcTimeNumerator(iOrderY + 1, iOrderX + 1) = -2*orderConstants(sumOrder + 1) * kernelConstants(iOrderX + 1) * kernelConstants(iOrderY + 1)/sampleCount;
            funcTimeExponent(iOrderY + 1, iOrderX + 1) = 1/(2 + sumOrder);
        end
    end
end
columnWeightedSpectra = cell(1, 6);
for iOrderX = 0:5
    columnWeightedSpectra{iOrderX + 1} = weightedSpectrum .* iPowers{iOrderX + 1};
end
baseOrderX = [5 4 3 2 1 0];
baseOrderY = [0 1 2 3 4 5];
objectiveCacheTimes = zeros(0, 1);
objectiveCacheValues = zeros(0, 1);
objectiveCacheMoments = zeros(0, 3);

tStar = smallestRoot(@objective, sampleCount);
[~, p02, p20, p11] = objectiveState(tStar);
tY = (p02^(3/4)/(4*pi*sampleCount*p20^(3/4)*(p11 + sqrt(p20*p02))))^(1/3);
tX = (p20^(3/4)/(4*pi*sampleCount*p02^(3/4)*(p11 + sqrt(p20*p02))))^(1/3);
bandwidth = sqrt([tX tY]).*scaling;

    function value = objective(t)
        value = objectiveState(t);
    end

    function [objectiveValue, p02, p20, p11] = objectiveState(t)
        iCached = find(objectiveCacheTimes == t, 1);
        if ~isempty(iCached)
            objectiveValue = objectiveCacheValues(iCached);
            p02 = objectiveCacheMoments(iCached, 1);
            p20 = objectiveCacheMoments(iCached, 2);
            p11 = objectiveCacheMoments(iCached, 3);
            return
        end

        [p02, p20, p11] = momentsForTime(t);
        sumFunc = p02 + p20 + 2 * p11;
        timeValue = (2*pi*sampleCount*sumFunc)^(-1/3);
        fixedPointValue = (t - timeValue)/timeValue;
        objectiveValue = t - fixedPointValue;
        objectiveCacheTimes(end + 1, 1) = t; %#ok<AGROW>
        objectiveCacheValues(end + 1, 1) = objectiveValue; %#ok<AGROW>
        objectiveCacheMoments(end + 1, :) = [p02 p20 p11]; %#ok<AGROW>
    end

    function [p02, p20, p11] = momentsForTime(t)
        baseValues = baseOrderFiveValues(t);
        f50 = baseValues(1);
        f41 = baseValues(2);
        f32 = baseValues(3);
        f23 = baseValues(4);
        f14 = baseValues(5);
        f05 = baseValues(6);

        f22 = funcFromChildren(2, 2, f32, f23);
        f13 = funcFromChildren(1, 3, f23, f14);
        f04 = funcFromChildren(0, 4, f14, f05);
        f40 = funcFromChildren(4, 0, f50, f41);
        f31 = funcFromChildren(3, 1, f41, f32);

        f12 = funcFromChildren(1, 2, f22, f13);
        f03 = funcFromChildren(0, 3, f13, f04);
        f30 = funcFromChildren(3, 0, f40, f31);
        f21 = funcFromChildren(2, 1, f31, f22);

        p02 = funcFromChildren(0, 2, f12, f03);
        p20 = funcFromChildren(2, 0, f30, f21);
        p11 = funcFromChildren(1, 1, f21, f12);
    end

    function value = funcFromChildren(orderX, orderY, leftChild, rightChild)
        timeValue = (funcTimeNumerator(orderY + 1, orderX + 1)/(leftChild + rightChild))^funcTimeExponent(orderY + 1, orderX + 1);
        value = psiValue(orderX, orderY, timeValue);
    end

    function values = baseOrderFiveValues(timeValue)
        weightedOrders = orderPowerMatrix .* exp(-iSquared * piSquared * timeValue);
        rowProducts = weightedOrders * weightedSpectrum;
        values = zeros(1, 6);
        for iValue = 1:6
            iOrderX = baseOrderX(iValue);
            iOrderY = baseOrderY(iValue);
            values(iValue) = psiScale(iOrderY + 1, iOrderX + 1) * ...
                (rowProducts(iOrderY + 1, :) * weightedOrders(iOrderX + 1, :).');
        end
    end

    function value = psiValue(orderX, orderY, timeValue)
        baseWeight = exp(-iSquared * piSquared * timeValue);
        weightedOrderY = iPowers{orderY + 1} .* baseWeight;
        value = psiScale(orderY + 1, orderX + 1) * ...
            (weightedOrderY * columnWeightedSpectra{orderX + 1} * baseWeight.');
    end
end

function value = kernelConstant(order)
value = (-1)^order * prod(1:2:(2*order - 1))/sqrt(2*pi);
end

function transformed = discreteCosineTransform1D(data)
nRows = size(data, 1);
weight = [1; 2 * exp(-1i * (1:(nRows - 1))' * pi/(2*nRows))];
transformed = real(weight .* fft([data(1:2:end, :); data(end:-2:2, :)]));
end

function transformed = discreteCosineTransform2D(data)
[nRows, nCols] = size(data);
if nRows ~= nCols
    error("KernelDensityEstimate:InvalidGridSize", ...
        "The internal two-dimensional Botev solver requires a square DCT grid.");
end

weight = [1; 2 * exp(-1i * (1:(nRows - 1))' * pi/(2*nRows))];
weight = weight(:, ones(1, nCols));
transformed = discreteCosineTransformColumns(discreteCosineTransformColumns(data).').';

    function output = discreteCosineTransformColumns(input)
        output = real(weight .* fft([input(1:2:end, :); input(end:-2:2, :)]));
    end
end

function histogramValues = normalizedHistogram2D(data, gridSize)
nRows = size(data, 1);
bins = zeros(nRows, size(data, 2));

for iDimension = 1:size(data, 2)
    [~, bins(:, iDimension)] = histc(data(:, iDimension), 0:1/gridSize:1, 1);
    bins(:, iDimension) = min(bins(:, iDimension), gridSize);
end

histogramValues = accumarray(bins(all(bins > 0, 2), :), 1/nRows, [gridSize gridSize]);
end

function value = smallestRoot(objective, sampleCount)
sampleCount = 50 * (sampleCount <= 50) + 1050 * (sampleCount >= 1050) + sampleCount * ((sampleCount < 1050) & (sampleCount > 50));
searchLimit = 1e-12 + 0.01 * (sampleCount - 50)/1000;
lowerBound = 0;
objectiveLower = objective(lowerBound);

while true
    objectiveUpper = objective(searchLimit);
    if hasRootBracket(objectiveLower, objectiveUpper)
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

function tf = hasRootBracket(lowerValue, upperValue)
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

function value = roundUpPowerOfTwo(value)
value = 2^ceil(log2(double(value)));
end
