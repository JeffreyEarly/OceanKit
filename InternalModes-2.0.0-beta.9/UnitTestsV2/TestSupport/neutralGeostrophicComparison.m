function metrics = neutralGeostrophicComparison(reference,candidate,z,weights)
% Compare normalized F/G pairs in positive physical-z H1 norms.
arguments
    reference struct
    candidate
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
truth = reference.evaluate(z);
if isstruct(candidate)
    values = candidate.evaluate(z);
    h = candidate.h;
    inverseDepth = candidate.inverseDepth;
    labels = candidate.modeNumber;
else
    values = sampleNeutralGeostrophicBasis(candidate,z);
    if isa(candidate,"IMGeostrophicZeroAPVModesBasis")
        h = reference.h;
        inverseDepth = reference.inverseDepth;
        labels = reference.modeNumber;
    else
        h = candidate.h;
        if isprop(candidate,"inverseDepth")
            inverseDepth = candidate.inverseDepth;
        else
            inverseDepth = 1./h;
        end
        labels = candidate.modeNumber;
    end
end
metrics.labelsMatch = isequal(labels,reference.modeNumber);
n = numel(reference.h);
if size(values.F,2) ~= n || ~metrics.labelsMatch
    metrics = failedComparison(metrics,n,"Mode counts or physical labels differ.");
    return
end
D = diff(reference.configuration.zDomain);
if reference.family ~= "zeroAPV"
    if reference.family == "mda"
        signs = sign(sum(weights.*truth.G.*values.G,1));
    else
        signs = sign(sum(weights.*truth.F.*values.F,1));
    end
    signs(signs == 0) = 1;
    for name = ["F","G","Fz","Gz"]
        values.(name) = values.(name).*signs;
    end
end
FNorm = h1Norm(truth.F,truth.Fz,weights,D);
GNorm = h1Norm(truth.G,truth.Gz,weights,D);
pairScale = max(FNorm,GNorm);
metrics.zeroFChannel = FNorm <= 1e3*eps.*max(pairScale,realmin);
metrics.zeroGChannel = GNorm <= 1e3*eps.*max(pairScale,realmin);
FNorm(metrics.zeroFChannel) = pairScale(metrics.zeroFChannel);
GNorm(metrics.zeroGChannel) = pairScale(metrics.zeroGChannel);
metrics.FH1Error = h1Norm(values.F-truth.F,values.Fz-truth.Fz,weights,D)./max(FNorm,realmin);
metrics.GH1Error = h1Norm(values.G-truth.G,values.Gz-truth.Gz,weights,D)./max(GNorm,realmin);
metrics.relativeHError = abs(h-reference.h)./max(abs(reference.h),realmin);
zero = isinf(reference.h);
metrics.relativeHError(zero) = double(~isinf(h(zero)) | inverseDepth(zero) ~= 0);
metrics.relativeInverseDepthError = abs(inverseDepth-reference.inverseDepth)./max(abs(reference.inverseDepth),realmin);
metrics.relativeInverseDepthError(zero) = double(inverseDepth(zero) ~= 0);
metrics.maximumFH1Error = max(metrics.FH1Error);
metrics.maximumGH1Error = max(metrics.GH1Error);
metrics.maximumRelativeHError = max(metrics.relativeHError);
metrics.maximumError = max([metrics.FH1Error metrics.GH1Error metrics.relativeHError]);
metrics.matchedModeCount = n;
metrics.reason = "";
metrics.sampledFGram = values.F.'*(weights.*values.F);
metrics.sampledGGram = values.G.'*(weights.*values.G);
end

function value = h1Norm(values,derivatives,weights,D)
value = sqrt(sum(weights.*(abs(values).^2+D^2*abs(derivatives).^2),1));
end

function metrics = failedComparison(metrics,n,reason)
metrics.zeroFChannel = false(1,n);
metrics.zeroGChannel = false(1,n);
metrics.FH1Error = Inf(1,n);
metrics.GH1Error = Inf(1,n);
metrics.relativeHError = Inf(1,n);
metrics.relativeInverseDepthError = Inf(1,n);
metrics.maximumFH1Error = Inf;
metrics.maximumGH1Error = Inf;
metrics.maximumRelativeHError = Inf;
metrics.maximumError = Inf;
metrics.matchedModeCount = 0;
metrics.reason = reason;
metrics.sampledFGram = [];
metrics.sampledGGram = [];
end
