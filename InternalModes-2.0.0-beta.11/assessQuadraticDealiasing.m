function report = assessQuadraticDealiasing(F,G,gridDegree,options)
% Measure modal bandwidth and apply an inexpensive quadratic-filtering policy.
%
% F and G contain one mode per column in the same supplied coordinate. With
% representation="values", rows are Chebyshev--Lobatto samples at
% $$\xi_q=-\cos(\pi q/(Q-1)),\quad q=0,\ldots,Q-1.$$
% With representation="coefficients", rows multiply Chebyshev polynomials
% of degree 0 through Q-1. A single row describes a constant. The caller
% evaluates every physical channel in the common coordinate, for example
% the simulation grid's affine WKB coordinate; coefficients in different
% native solver coordinates must not be compared directly. gridDegree is
% the simulation polynomial degree, even when the supplied samples are finer.
%
% The effective degree is the smallest d containing energyFraction of the
% Chebyshev-weighted energy, proportional to
% $$2|a_0|^2+\sum_{n=1}^{Q-1}|a_n|^2.$$
% Full energy coverage includes every nonzero coefficient, even when its
% squared amplitude is below floating-point energy resolution. Channels are
% scaled independently before squaring. An absent channel [] or
% an exactly zero column contributes degree and tail fraction zero. Each
% mode's degree is the larger of its F and G degrees. Effective-bandwidth
% acceptance requires this degree <= floor(bandwidthFraction*gridDegree).
% This is a shape heuristic, not a bound on nonlinear aliasing or a check
% that the supplied sample resolution is sufficient.
%
% The fixed-fraction policy accepts the first floor(retainedFraction*M)
% columns. M must describe the full linearly accepted prefix, before any
% explicit requested count is applied. The none policy accepts all columns.
% Both policies skip bandwidth measurement and return NaN spectral metrics.
% Per-mode masks may contain gaps under effectiveBandwidth; the caller owns
% mode ordering, contiguous-prefix selection, and explicit-count policy.
%
% Values cost O(M Q log Q), coefficients cost O(M Q), and additional working
% storage is O(Q min(M,batchSize)) plus O(M) reporting. Caller-owned samples
% are excluded. No eigensolves or mode-pair products are performed.
%
% - Topic: Analyze modes
% - Declaration: report = assessQuadraticDealiasing(F,G,gridDegree,options)
% - Parameter F: real Q-by-M common-coordinate values or coefficients; [] omits F
% - Parameter G: real Q-by-M common-coordinate values or coefficients; [] omits G
% - Parameter gridDegree: nonnegative simulation polynomial degree
% - Parameter options.quadraticDealiasing: "none", "fixedFraction", or "effectiveBandwidth"
% - Parameter options.retainedFraction: fixed-fraction prefix fraction in [0,1]
% - Parameter options.energyFraction: effective-bandwidth energy coverage in (0,1]
% - Parameter options.bandwidthFraction: fraction of gridDegree allowed in [0,1]
% - Parameter options.representation: "values" or "coefficients"
% - Parameter options.batchSize: maximum number of columns transformed together
% - Returns report: policy, limits, column-oriented masks, effective degrees, and tail energy fractions
arguments (Input)
    F (:,:) double {mustBeReal,mustBeFinite}
    G (:,:) double {mustBeReal,mustBeFinite}
    gridDegree (1,1) double {mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
    options.quadraticDealiasing (1,1) string {mustBeMember(options.quadraticDealiasing,["none","fixedFraction","effectiveBandwidth"])} = "effectiveBandwidth"
    options.retainedFraction (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative,mustBeLessThanOrEqual(options.retainedFraction,1)} = 2/3
    options.energyFraction (1,1) double {mustBeReal,mustBeFinite,mustBePositive,mustBeLessThanOrEqual(options.energyFraction,1)} = .99
    options.bandwidthFraction (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative,mustBeLessThanOrEqual(options.bandwidthFraction,1)} = 2/3
    options.representation (1,1) string {mustBeMember(options.representation,["values","coefficients"])} = "values"
    options.batchSize (1,1) double {mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 128
end
arguments (Output)
    report (1,1) struct
end

hasF = ~isequal(size(F),[0 0]);
hasG = ~isequal(size(G),[0 0]);
if hasF && hasG && ~isequal(size(F),size(G))
    error("assessQuadraticDealiasing:ChannelSizeMismatch","F and G must have the same sample and mode dimensions; use [] to omit a channel.");
end
if hasF
    [sampleCount,nModes] = size(F);
else
    [sampleCount,nModes] = size(G);
end
if sampleCount == 0 && nModes > 0
    error("assessQuadraticDealiasing:MissingSamples","A nonempty mode inventory requires at least one sample or coefficient per column.");
end

report = struct(quadraticDealiasing=options.quadraticDealiasing,gridDegree=gridDegree,nModes=nModes,sampleCount=sampleCount,retainedFraction=options.retainedFraction,energyFraction=options.energyFraction,bandwidthFraction=options.bandwidthFraction,bandwidthLimit=floor(options.bandwidthFraction*gridDegree),representation=options.representation);
report.accepted = true(nModes,1);
report.effectiveDegreeF = nan(nModes,1);
report.effectiveDegreeG = nan(nModes,1);
report.effectiveDegree = nan(nModes,1);
report.tailEnergyFractionF = nan(nModes,1);
report.tailEnergyFractionG = nan(nModes,1);
if options.quadraticDealiasing == "fixedFraction"
    report.accepted(floor(options.retainedFraction*nModes)+1:end) = false;
elseif options.quadraticDealiasing == "effectiveBandwidth"
    [report.effectiveDegreeF,report.tailEnergyFractionF] = measureChannel(F,nModes,report.bandwidthLimit,options);
    [report.effectiveDegreeG,report.tailEnergyFractionG] = measureChannel(G,nModes,report.bandwidthLimit,options);
    report.effectiveDegree = max(report.effectiveDegreeF,report.effectiveDegreeG);
    report.accepted = report.effectiveDegree <= report.bandwidthLimit;
end
end

function [degree,tailFraction] = measureChannel(samples,nModes,degreeLimit,options)
degree = zeros(nModes,1);
tailFraction = zeros(nModes,1);
if isempty(samples)
    return
end
sampleCount = size(samples,1);
for first = 1:options.batchSize:nModes
    columns = first:min(first+options.batchSize-1,nModes);
    coefficients = samples(:,columns);
    scale = max(abs(coefficients),[],1);
    scale(scale == 0) = 1;
    coefficients = coefficients./scale;
    if options.representation == "values" && sampleCount > 1
        % Even extension gives the DCT-I using only base MATLAB's FFT.
        coefficients = flipud(coefficients);
        coefficients = real(fft([coefficients;coefficients(end-1:-1:2,:)],[],1))/(sampleCount-1);
        coefficients = coefficients(1:sampleCount,:);
        coefficients([1 end],:) = coefficients([1 end],:)/2;
    end
    energy = coefficients.^2;
    energy(1,:) = 2*energy(1,:);
    cumulative = cumsum(energy,1);
    total = cumulative(end,:);
    if options.energyFraction == 1
        % Inspect coefficients before squaring; tiny tails may disappear
        % from accumulated energy or underflow when squared.
        if options.representation == "coefficients"
            nonzero = samples(:,columns) ~= 0;
        else
            nonzero = coefficients ~= 0;
        end
        degree(columns) = max((0:sampleCount-1)'.*nonzero,[],1).';
    else
        degree(columns) = sum(cumulative < options.energyFraction*total,1).';
    end
    denominator = total;
    denominator(denominator == 0) = 1;
    tailFraction(columns) = (sum(energy(degreeLimit+2:end,:),1)./denominator).';
end
end
