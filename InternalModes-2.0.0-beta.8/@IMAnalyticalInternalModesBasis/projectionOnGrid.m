function projection = projectionOnGrid(self,z,weights,options)
% Prepare physical scalar pairings on the supplied fixed rule.
%
% Select distinct columns before solving the retained physical form.
% Omitted columns select all; [] selects none. Physical preparation does not
% infer missing independent endpoint or companion observations.
%
% The physical dual applies the continuous coefficient normalization to
% sampled pairings. Sampled-overlap recovery is a separate operation.
%
% - Topic: Build discrete transforms
% - Parameter z: increasing physical sample column
% - Parameter weights: fixed physical quadrature weights
% - Parameter options.variable: required independent scalar channel
% - Parameter options.columns: distinct column indices, all by default
% - Parameter options.prepareRecovery: also prepare sampled recovery, false by default
% - Parameter options.failureMode: diagnostic or error for unresolved recovery
% - Parameter options.rankTolerance: relative numerical rank cutoff
% - Returns projection: physical projection record
arguments (Input)
    self (1,1) IMAnalyticalInternalModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:self.nColumns
    options.prepareRecovery (1,1) logical = false
    options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
    options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
end
arguments (Output)
    projection (1,1) IMProjection
end
args = namedargs2cell(options);
projection = IMScalarProjectionTools.projectionOnGrid(self,z,weights,args{:});
end
