function projection = fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,options)
% Construct a prescribed physical source dual without inventing synthesis data.
%
% `samplePairingMatrix` maps sampled source values into signed coefficient
% pairings. `sampleGram` supplies the signed coefficient system to solve;
% it need not equal a single observed component's Gram. `targetGram` is the
% corresponding continuous reference system. Prefixes restrict rows and
% coefficient systems together and solve each restricted system separately.
%
% The caller owns the physical recipe and supplies explicit provenance.
% The object can measure product errors using qualified reference pairings,
% but Gram and round-trip assessment are unavailable without a synthesis
% basis. Their diagnostics are NaN and capability flags are false.
%
% - Topic: Create projections
% - Declaration: projection = IMProjection.fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,options)
% - Parameter samplePairingMatrix: finite real or complex column-by-sample pairing operator
% - Parameter sampleGram: real symmetric signed system for sampled coefficients
% - Parameter targetGram: real symmetric continuous coefficient system
% - Parameter options.majorantGramMatrix: positive coefficient error metric
% - Parameter options.activeColumnMask: active output coordinates, default all
% - Parameter options.columnLabels: distinct scientific output labels
% - Parameter options.provenance: required description of the supplied physical dual
% - Returns projection: source projection with no sampled synthesis basis
arguments (Input)
    samplePairingMatrix (:,:) double {mustBeFinite}
    sampleGram (:,:) double {mustBeReal,mustBeFinite}
    targetGram (:,:) double {mustBeReal,mustBeFinite}
    options.majorantGramMatrix (:,:) double {mustBeReal,mustBeFinite} = []
    options.activeColumnMask (1,:) logical = true(1,size(samplePairingMatrix,1))
    options.columnLabels (1,:) string = string(1:size(samplePairingMatrix,1))
    options.provenance (1,1) struct
end
arguments (Output)
    projection (1,1) IMProjection
end
if isempty(fieldnames(options.provenance))
    error("IMProjection:MissingProvenance","A prescribed physical dual requires explicit nonempty provenance.");
end
projection = IMProjection();
projection.projectionKind = "prescribedDual";
projection = projection.initializeProjection(samplePairingMatrix,sampleGram,targetGram,options.majorantGramMatrix,options.activeColumnMask,options.columnLabels);
projection.projectionKind = "prescribedDual";
projection.supportsGramAssessment = false;
projection.supportsRoundTripAssessment = false;
projection.provenance = options.provenance;
projection.sampledBasis = [];
projection.metricMatrix = [];
projection.gramError = NaN;
projection.roundTripError = NaN;
projection.inverseMatrixConditionNumber = NaN;
end
