function result = prepareIMRecovery(A,options)
% Prepare strict sampled recovery without choosing modes, grids or quadrature.
%
% On the declared active columns, C maps dimensionless coefficients to the
% supplied coefficient coordinates and D whitens the positive observation
% norm. Least squares uses the independently declared fitting factor F:
% $$R=C(FAC)^{+}F.$$
% The full-column-rank solve uses the weighted synthesis directly, never its
% normal matrix. Corrected physical recovery retains the supplied dual B and
% solves in an equilibrated orthonormal frame of
% $$B_s=C^{-1}BD^{-1},\qquad A_s=DAC.$$
% Its result is the same mathematical map as $$(BA)^{-1}B$$. A missing D is
% permitted for algebraic adapters, but observation/noise/geometric evidence
% remains unavailable. F-only least squares likewise declares no noise norm.
% Factorizations use bounded matrices and binary scale factors. Spectra that
% exceed floating-point range retain scaled singular values and their scale
% factors; rank never depends on overflowing or underflowing singular values.
% result.specification retains the resolved declared factors, including their
% absence, independently of internal numerical preconditioning.
%
% Inactive coordinates must have identically zero synthesis columns; their
% aligned forward rows remain exactly zero. Failure in diagnostic mode returns
% no applicable forward matrix. Neither policy silently truncates directions.
%
% - Topic: Developer topics — Prepare recovery
% - Declaration: result = prepareIMRecovery(A,options)
% - Parameter A: finite complex sample-by-column synthesis matrix
% - Parameter options.policy: sampledLeastSquares, correctedPhysicalDual or extendedSpan
% - Parameter options.columnLabels: distinct nonempty scientific string labels
% - Parameter options.activeColumnMask: declared recoverable coefficient columns
% - Parameter options.coefficientScale: nonsingular active C; a vector denotes its diagonal
% - Parameter options.observationFactor: nonsingular sample-space D; a vector denotes its diagonal
% - Parameter options.fitFactor: rectangular F; defaults to D when supplied
% - Parameter options.physicalDual: unchanged full-aligned coefficient-by-sample B
% - Parameter options.rankTolerance: relative singular-value cutoff, strictly between zero and one
% - Parameter options.failureMode: error or diagnostic
% - Parameter options.extendedSynthesis: larger sample-by-column synthesis for extendedSpan
% - Parameter options.extendedColumnLabels: scientific labels of the larger span
% - Parameter options.extendedActiveColumnMask: activity of the larger span
% - Parameter options.extendedCoefficientScale: C on the larger active span
% - Parameter options.retainedColumnIndices: larger-span indices matching A exactly
% - Parameter options.extendedPolicy: recovery policy on the larger span
% - Parameter options.extendedPhysicalDual: physical B aligned to the larger span
% - Returns result: recovery matrix, unchanged physical dual and separate numerical evidence
% - Developer: true
arguments (Input)
    A (:,:) double {mustBeFinite}
    options.policy (1,1) string {mustBeMember(options.policy,["sampledLeastSquares","correctedPhysicalDual","extendedSpan"])} = "sampledLeastSquares"
    options.columnLabels (1,:) string = string(1:size(A,2))
    options.activeColumnMask (1,:) logical = true(1,size(A,2))
    options.coefficientScale (:,:) double {mustBeFinite} = []
    options.observationFactor (:,:) double {mustBeFinite} = []
    options.fitFactor (:,:) double {mustBeFinite} = []
    options.physicalDual (:,:) double {mustBeFinite} = []
    options.rankTolerance (1,1) double {mustBeReal,mustBePositive,mustBeFinite} = 1e-12
    options.failureMode (1,1) string {mustBeMember(options.failureMode,["error","diagnostic"])} = "error"
    options.extendedSynthesis (:,:) double {mustBeFinite} = []
    options.extendedColumnLabels (1,:) string = strings(1,0)
    options.extendedActiveColumnMask (1,:) logical = false(1,0)
    options.extendedCoefficientScale (:,:) double {mustBeFinite} = []
    options.retainedColumnIndices (1,:) double {mustBeInteger,mustBePositive,mustBeFinite} = zeros(1,0)
    options.extendedPolicy (1,1) string {mustBeMember(options.extendedPolicy,["sampledLeastSquares","correctedPhysicalDual"])} = "sampledLeastSquares"
    options.extendedPhysicalDual (:,:) double {mustBeFinite} = []
end
arguments (Output)
    result (1,1) struct
end
[m,n] = size(A);
if options.rankTolerance >= 1
    error("IMRecovery:InvalidRankTolerance","Use a relative rank tolerance strictly between zero and one.");
end
if numel(options.columnLabels) ~= n || any(ismissing(options.columnLabels) | strlength(options.columnLabels) == 0) || numel(unique(options.columnLabels)) ~= n
    error("IMRecovery:InvalidColumnLabels","Supply one distinct nonempty scientific label per synthesis column.");
end
if numel(options.activeColumnMask) ~= n
    error("IMRecovery:InvalidActivity","Supply one declared activity flag per synthesis column.");
end
if any(A(:,~options.activeColumnMask) ~= 0,"all")
    error("IMRecovery:InvalidInactiveColumn","Declared inactive synthesis columns must be identically zero; sampled-zero active columns remain active.");
end
if ~isempty(options.physicalDual) && ~isequal(size(options.physicalDual),[n m])
    error("IMRecovery:InvalidPhysicalDual","Supply physicalDual in the full aligned coefficient-by-sample layout.");
end
active = find(options.activeColumnMask);
na = numel(active);
C = coordinateFactor(options.coefficientScale,na,"coefficientScale");
hasObservationNorm = ~isempty(options.observationFactor);
if hasObservationNorm
    D = coordinateFactor(options.observationFactor,m,"observationFactor");
else
    D = [];
end
F = options.fitFactor;
hasFitFactor = ~isequal(size(F),[0 0]);
if ~hasFitFactor && hasObservationNorm
    F = D;
    hasFitFactor = true;
end
if hasFitFactor && size(F,2) ~= m
    error("IMRecovery:InvalidFitFactor","fitFactor must have one column per observation; zero fitting directions are permitted.");
end
[An,aMagnitude] = boundedMatrix(A(:,active));
[Cn,cMagnitude] = boundedMatrix(C);
if hasObservationNorm
    [Dn,dMagnitude] = boundedMatrix(D);
    [As,asMagnitude] = boundedMatrix(Dn*An*Cn);
else
    Dn = [];
    dMagnitude = 1;
    [As,asMagnitude] = boundedMatrix(An*Cn);
end
unavailable = unavailableEvidence(options.rankTolerance);
diagnostics = struct(observation=unavailable,fit=unavailable,test=unavailable,overlap=unavailable,roundTripError=NaN,scaledRoundTripError=NaN,roundTripStatus="unavailable",definingEquationResidual=NaN,definingEquationResidualStatus="unavailable",scaledDefiningEquationResidual=NaN,observationNoiseGain=NaN,observationNoiseGainStatus="unavailable",noiseOptimalForObservationNorm=false,extendedRejectionError=NaN,principalCosines=zeros(0,1),supportsPrincipalGeometry=false,extended=struct);
if hasObservationNorm, diagnostics.observation = matrixEvidence(As,na,options.rankTolerance,[dMagnitude aMagnitude cMagnitude asMagnitude]); end
declaredD = [];
if hasObservationNorm, declaredD = D; end
specification = struct(policy=options.policy,columnLabels=options.columnLabels,activeColumnMask=options.activeColumnMask,coefficientScale=C,observationFactor=declaredD,fitFactor=F,rankTolerance=options.rankTolerance,failureMode=options.failureMode,extendedSynthesis=options.extendedSynthesis,extendedColumnLabels=options.extendedColumnLabels,extendedActiveColumnMask=options.extendedActiveColumnMask,extendedCoefficientScale=options.extendedCoefficientScale,retainedColumnIndices=options.retainedColumnIndices,extendedPolicy=options.extendedPolicy,extendedPhysicalDual=options.extendedPhysicalDual);
result = struct(forwardMatrix=[],physicalDual=options.physicalDual,success=false,status="unprepared",policy=options.policy,columnLabels=options.columnLabels,activeColumnMask=options.activeColumnMask,diagnostics=diagnostics,specification=specification);
if options.policy == "extendedSpan"
    result = extendedRecovery(A,options,C,D,hasObservationNorm,result);
    return
end
if na == 0
    result.forwardMatrix = zeros(n,m);
    result.success = true;
    result.status = "emptyActiveSpan";
    result.diagnostics.roundTripError = 0;
    result.diagnostics.scaledRoundTripError = 0;
    result.diagnostics.roundTripStatus = "measured";
    if hasObservationNorm
        result.diagnostics.observationNoiseGain = 0;
        result.diagnostics.observationNoiseGainStatus = "measured";
    end
    if options.policy == "sampledLeastSquares" && hasFitFactor
        result.diagnostics.fit = matrixEvidence(zeros(size(F,1),0),0,options.rankTolerance);
        result.diagnostics.noiseOptimalForObservationNorm = hasObservationNorm && isequal(F,D);
    elseif options.policy == "correctedPhysicalDual"
        result.diagnostics.test = matrixEvidence(zeros(m,0),0,options.rankTolerance);
        result.diagnostics.overlap = matrixEvidence(zeros(0),0,options.rankTolerance);
        result.diagnostics.definingEquationResidual = 0;
        result.diagnostics.definingEquationResidualStatus = "measured";
        result.diagnostics.scaledDefiningEquationResidual = 0;
    end
    return
end
if hasObservationNorm && result.diagnostics.observation.rank < na
    result = rejectRecovery(result,"sampleRankDeficient",options.failureMode,"The declared observation norm does not resolve every active column at the stated rank tolerance.");
    return
end
if options.policy == "sampledLeastSquares"
    if ~hasFitFactor
        error("IMRecovery:MissingFitFactor","Declare observationFactor or fitFactor for sampled least squares.");
    end
    [Fn,fMagnitude] = boundedMatrix(F);
    [fitted,fitMagnitude] = boundedMatrix(Fn*An*Cn);
    result.diagnostics.fit = matrixEvidence(fitted,na,options.rankTolerance,[fMagnitude aMagnitude cMagnitude fitMagnitude]);
    if result.diagnostics.fit.rank < na
        result = rejectRecovery(result,"fitRankDeficient",options.failureMode,"The fitting factor does not resolve every active column; provide more fitting information or request a different explicit span.");
        return
    end
    [Q,R] = qr(fitted,0);
    forward = scaleByFactors(Cn*(R\(Q'*Fn)),[],[aMagnitude fitMagnitude]);
    result.diagnostics.noiseOptimalForObservationNorm = hasObservationNorm && isequal(F,D);
else
    if isempty(options.physicalDual)
        error("IMRecovery:MissingPhysicalDual","Supply the full aligned physicalDual for corrected physical recovery.");
    end
    B = options.physicalDual(active,:);
    % C^{-1} is an invertible left change and does not change the test row
    % space. Normalize physical rows first, before any potentially enormous
    % row norm or global coefficient scaling can lose their directions.
    [boundedRows,rowMagnitudes,rowLengths] = boundedRowsOf(B);
    whitenedRows = boundedRows;
    if hasObservationNorm, whitenedRows = boundedRows/Dn; end
    if any(~isfinite(whitenedRows),"all")
        result = rejectRecovery(result,"unrepresentableIntermediate",options.failureMode,"Whitening physical test rows exceeds floating-point range in these coordinates; revise the declared coordinates or use a separately requested recovery policy.");
        return
    end
    [equilibrated,~,~] = boundedRowsOf(whitenedRows);
    result.diagnostics.test = matrixEvidence(equilibrated',na,options.rankTolerance);
    result.diagnostics.test.rowNorms = rowMagnitudes.*rowLengths;
    result.diagnostics.test.rowMagnitudeScales = rowMagnitudes;
    result.diagnostics.test.rowNormalizedLengths = rowLengths;
    result.diagnostics.test.rowNormsCoordinateSystem = "physicalDualBeforeWhitening";
    if result.diagnostics.test.rank < na
        result = rejectRecovery(result,"testRankDeficient",options.failureMode,"The equilibrated physical test space lacks an independent row for every active coefficient.");
        return
    end
    [testFrame,~,~] = svd(equilibrated',"econ");
    rowFrame = testFrame(:,1:na)';
    overlap = rowFrame*As;
    result.diagnostics.overlap = matrixEvidence(overlap,na,options.rankTolerance,[dMagnitude aMagnitude cMagnitude asMagnitude]);
    [Bn,bMagnitude] = boundedMatrix(B);
    result.diagnostics.overlap.raw = matrixEvidence((Cn\(Bn*An))*Cn,na,options.rankTolerance,[bMagnitude aMagnitude]);
    if hasObservationNorm
        [trialFrame,~] = qr(As,0);
        result.diagnostics.principalCosines = svd(rowFrame*trialFrame);
        result.diagnostics.supportsPrincipalGeometry = true;
    end
    if result.diagnostics.overlap.rank < na
        result = rejectRecovery(result,"overlapRankDeficient",options.failureMode,"The physical trial/test overlap is deficient after row-space preparation; choose a different explicit span or policy.");
        return
    end
    scaledForward = overlap\rowFrame;
    forwardBody = Cn*scaledForward;
    if hasObservationNorm, forwardBody = forwardBody*Dn; end
    forward = scaleByFactors(forwardBody,[],[aMagnitude asMagnitude]);
    result.diagnostics.scaledDefiningEquationResidual = norm((equilibrated*As)*scaledForward-equilibrated,2)/max(norm(equilibrated,2),realmin);
    if all(isfinite(forward),"all")
        [result.diagnostics.definingEquationResidual,result.diagnostics.definingEquationResidualStatus] = weakDefiningResidual(B,A(:,active),forward);
    end
end
result.forwardMatrix = zeros(n,m,"like",forward);
result.forwardMatrix(active,:) = forward;
result = finishRecovery(result,A,C,D,active,hasObservationNorm);
end

function result = extendedRecovery(A,options,C,D,hasObservationNorm,result)
[m,n] = size(A);
extendedA = options.extendedSynthesis;
ne = size(extendedA,2);
indices = options.retainedColumnIndices;
if size(extendedA,1) ~= m || numel(indices) ~= n || numel(unique(indices)) ~= n || any(indices > ne)
    error("IMRecovery:InvalidExtendedSpan","Supply an extended synthesis on the same observations and one distinct valid retained index per synthesis column.");
end
labels = options.extendedColumnLabels;
mask = options.extendedActiveColumnMask;
if numel(labels) ~= ne || numel(mask) ~= ne || ~isequal(labels(indices),options.columnLabels) || ~isequal(mask(indices),options.activeColumnMask) || ~isequal(extendedA(:,indices),A)
    error("IMRecovery:InconsistentExtendedSpan","Retained columns, labels and activity must exactly match their declared extended-span entries.");
end
prepared = prepareIMRecovery(extendedA,policy=options.extendedPolicy,columnLabels=labels,activeColumnMask=mask,coefficientScale=options.extendedCoefficientScale,observationFactor=options.observationFactor,fitFactor=options.fitFactor,physicalDual=options.extendedPhysicalDual,rankTolerance=options.rankTolerance,failureMode="diagnostic");
result.diagnostics.extended = prepared.diagnostics;
result.specification.extendedCoefficientScale = prepared.specification.coefficientScale;
result.specification.extendedRecoverySpecification = prepared.specification;
if ~prepared.success
    result = rejectRecovery(result,"extendedRecoveryFailed",options.failureMode,"The declared extended span failed its requested recovery policy; inspect diagnostics.extended.");
    return
end
result.forwardMatrix = prepared.forwardMatrix(indices,:);
omitted = mask & ~ismember(1:ne,indices);
result.diagnostics.extendedRejectionError = stableProductNorm(result.forwardMatrix,extendedA(:,omitted));
result = finishRecovery(result,A,C,D,options.activeColumnMask,hasObservationNorm);
end

function result = finishRecovery(result,A,C,D,active,hasObservationNorm)
n = size(A,2);
if any(~isfinite(result.forwardMatrix),"all")
    result = rejectRecovery(result,"unrepresentableRecovery",result.specification.failureMode,"The requested recovery operator is outside the floating-point range; revise the declared coordinates or observations.");
    return
end
if any(all(result.forwardMatrix(active,:) == 0,2))
    result = rejectRecovery(result,"unrepresentableRecovery",result.specification.failureMode,"An active recovery row underflowed to zero; revise the declared coordinates or observations.");
    return
end
delta = stableProduct(result.forwardMatrix,A)-diag(double(result.activeColumnMask));
result.diagnostics.roundTripError = stableNorm(delta);
[Cn,~] = boundedMatrix(C);
result.diagnostics.scaledRoundTripError = stableNorm((Cn\delta(active,active))*Cn);
result.diagnostics.roundTripStatus = rangeStatus(max(result.diagnostics.roundTripError,result.diagnostics.scaledRoundTripError));
if hasObservationNorm
    [Rn,rMagnitude] = boundedMatrix(result.forwardMatrix(active,:));
    [Cn,cMagnitude] = boundedMatrix(C);
    [Dn,dMagnitude] = boundedMatrix(D);
    noiseBody = (Cn\Rn)/Dn;
    noiseNorm = stableNorm(noiseBody);
    result.diagnostics.observationNoiseGain = scaleByFactors(noiseNorm,rMagnitude,[cMagnitude dMagnitude]);
    result.diagnostics.observationNoiseGainStatus = rangeStatus(result.diagnostics.observationNoiseGain);
    if noiseNorm > 0 && result.diagnostics.observationNoiseGain == 0
        result.diagnostics.observationNoiseGainStatus = "outsideFloatingPointRange";
    end
end
result.success = true;
result.status = "ready";
if n == 0, result.status = "emptyActiveSpan"; end
end

function factor = coordinateFactor(value,count,name)
if isempty(value)
    factor = eye(count);
elseif isvector(value) && numel(value) == count
    factor = diag(value);
else
    factor = value;
end
if ~isequal(size(factor),[count count])
    error("IMRecovery:InvalidCoordinateFactor","%s must be nonsingular on its declared coordinates; coefficientScale uses active coordinates only.",name);
end
if count > 0 && rcond(factor) == 0
    error("IMRecovery:SingularCoordinateFactor","%s must be nonsingular; zero fitting directions belong in fitFactor.",name);
end
end

function evidence = matrixEvidence(matrix,columnCount,relativeTolerance,scaleFactors)
if nargin < 4, scaleFactors = 1; end
if any(~isfinite(matrix),"all")
    evidence = unavailableEvidence(relativeTolerance);
    evidence.status = "outsideFloatingPointRange";
    return
end
[unitMatrix,magnitude] = boundedMatrix(matrix);
scaledSingularValues = svd(unitMatrix,"econ");
scaleFactors = [scaleFactors magnitude];
singularValues = scaleByFactors(scaledSingularValues,scaleFactors,[]);
threshold = 0;
if ~isempty(scaledSingularValues), threshold = relativeTolerance*scaledSingularValues(1); end
rankValue = nnz(scaledSingularValues > threshold);
conditionNumber = Inf;
if columnCount == 0
    conditionNumber = 1;
elseif rankValue == columnCount
    conditionNumber = scaledSingularValues(1)/scaledSingularValues(columnCount);
end
representable = all(isfinite(singularValues)) && all(singularValues(scaledSingularValues > 0) > 0);
status = "measured";
if ~representable, status = "scaledSpectrum"; end
evidence = struct(available=true,status=status,rank=rankValue,columnCount=columnCount,conditionNumber=conditionNumber,singularValues=singularValues,scaledSingularValues=scaledSingularValues,singularValueScaleFactors=scaleFactors,singularValuesRepresentable=representable,relativeTolerance=relativeTolerance,absoluteTolerance=scaleByFactors(threshold,scaleFactors,[]));
end

function evidence = unavailableEvidence(relativeTolerance)
evidence = struct(available=false,status="unavailable",rank=NaN,columnCount=NaN,conditionNumber=NaN,singularValues=zeros(0,1),scaledSingularValues=zeros(0,1),singularValueScaleFactors=zeros(1,0),singularValuesRepresentable=false,relativeTolerance=relativeTolerance,absoluteTolerance=NaN);
end

function [unitMatrix,magnitude] = boundedMatrix(matrix)
% Component maxima remain representable even when a complex modulus or a
% Euclidean norm exceeds realmax. Dividing before summation protects both.
if any(~isfinite(matrix),"all")
    error("IMRecovery:UnrepresentableIntermediate","A scaled intermediate is nonfinite; revise the declared coordinate factors.");
end
magnitude = max([abs(real(matrix(:)));abs(imag(matrix(:)));0]);
if magnitude == 0, magnitude = 1; end
unitMatrix = matrix./magnitude;
end

function [unitRows,magnitudes,lengths] = boundedRowsOf(matrix)
if any(~isfinite(matrix),"all")
    error("IMRecovery:UnrepresentableIntermediate","Whitening the physical test rows produced a nonfinite intermediate.");
end
magnitudes = max(max(abs(real(matrix)),abs(imag(matrix))),[],2);
divisors = magnitudes;
divisors(divisors == 0) = 1;
unitRows = matrix./divisors;
lengths = vecnorm(unitRows,2,2);
divisors = lengths;
divisors(divisors == 0) = 1;
unitRows = unitRows./divisors;
end

function value = scaleByFactors(value,numeratorFactors,denominatorFactors)
% Binary exponents avoid first forming an overflowing or underflowing scale
% whose product with the requested output is nevertheless representable.
[numerators,numeratorExponents] = log2(numeratorFactors);
[denominators,denominatorExponents] = log2(denominatorFactors);
fraction = prod(numerators)/prod(denominators);
[fraction,exponent] = log2(fraction);
exponent = exponent+sum(numeratorExponents)-sum(denominatorExponents);
if isreal(value)
    value = scaleRealComponent(value,fraction,exponent);
else
    value = complex(scaleRealComponent(real(value),fraction,exponent),scaleRealComponent(imag(value),fraction,exponent));
end
end

function value = scaleRealComponent(value,fraction,exponent)
[valueFraction,valueExponent] = log2(value);
value = pow2(valueFraction.*fraction,valueExponent+exponent);
end

function product = stableProduct(left,right)
[unitLeft,leftMagnitude] = boundedMatrix(left);
[unitRight,rightMagnitude] = boundedMatrix(right);
product = scaleByFactors(unitLeft*unitRight,[leftMagnitude rightMagnitude],[]);
end

function [value,status] = stableProductNorm(left,right)
if any(~isfinite(left),"all") || any(~isfinite(right),"all")
    value = Inf;
    status = "outsideFloatingPointRange";
    return
end
[unitLeft,leftMagnitude] = boundedMatrix(left);
[unitRight,rightMagnitude] = boundedMatrix(right);
unitNorm = stableNorm(unitLeft*unitRight);
value = scaleByFactors(unitNorm,[leftMagnitude rightMagnitude],[]);
status = rangeStatus(value);
if unitNorm > 0 && value == 0, status = "outsideFloatingPointRange"; end
end

function [value,status] = weakDefiningResidual(B,A,R)
% Associate as (B*A)*R-B: observations may be large, while the coefficient
% dimension is small. Never form the observation-space projector A*R.
[Bn,bMagnitude] = boundedMatrix(B);
[An,aMagnitude] = boundedMatrix(A);
[Rn,rMagnitude] = boundedMatrix(R);
residual = scaleByFactors((Bn*An)*Rn,[aMagnitude rMagnitude],[])-Bn;
unitNorm = stableNorm(residual);
value = scaleByFactors(unitNorm,bMagnitude,[]);
status = rangeStatus(value);
if unitNorm > 0 && value == 0, status = "outsideFloatingPointRange"; end
end

function value = stableNorm(matrix)
if any(~isfinite(matrix),"all")
    value = Inf;
    return
end
[unitMatrix,magnitude] = boundedMatrix(matrix);
value = scaleByFactors(norm(unitMatrix,2),magnitude,[]);
end

function status = rangeStatus(value)
status = "measured";
if ~isfinite(value), status = "outsideFloatingPointRange"; end
end

function result = rejectRecovery(result,status,failureMode,message)
result.success = false;
result.status = status;
result.forwardMatrix = [];
if failureMode == "error"
    identifier = "IMRecovery:"+upper(extractBefore(status,2))+extractAfter(status,1);
    error(identifier,"%s",message);
end
end
