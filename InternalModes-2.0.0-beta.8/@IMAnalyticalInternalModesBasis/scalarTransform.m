function transform = scalarTransform(self,z,weights,options)
% Prepare scalar recovery separately from physical source projection.
%
% The physical source dual remains available on transform.projection.
% Numerical recovery uses an explicit policy and never substitutes a
% rank-deficient inverse. No new solve, point rule or quadrature fit occurs.
%
% The default corrected-physical-dual map is $$R=(BA)^{-1}B$$ for scalar
% synthesis $$A$$ and continuous-normalized physical dual $$B$$. The
% uncorrected $$B$$ remains available independently of successful recovery.
%
% - Topic: Build discrete transforms
% - Parameter z: increasing physical sample column
% - Parameter weights: fixed physical quadrature weights
% - Parameter options.variable: required independent scalar channel
% - Parameter options.recovery: optional aligned prepared IMDiscreteTransform
% - Parameter options.policy: correctedPhysicalDual or sampledLeastSquares
% - Parameter options.coefficientScale: nonsingular active amplitude scaling
% - Parameter options.observationFactor: declared positive observation whitening
% - Parameter options.fitFactor: explicit least-squares fitting factor
% - Parameter options.rankTolerance: relative numerical rank cutoff
% - Parameter options.failureMode: diagnostic or error for unresolved recovery
% - Returns transform: physical projection record and independently prepared recovery
arguments (Input)
    self (1,1) IMAnalyticalInternalModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.recovery = []
    options.policy (1,1) string {mustBeMember(options.policy,["correctedPhysicalDual","sampledLeastSquares"])} = "correctedPhysicalDual"
    options.coefficientScale (:,:) double {mustBeFinite} = []
    options.observationFactor (:,:) double {mustBeFinite} = []
    options.fitFactor (:,:) double {mustBeFinite} = []
    options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
    options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
end
arguments (Output)
    transform (1,1) IMScalarDiscreteTransform
end
args = namedargs2cell(options);
transform = IMScalarProjectionTools.scalarTransform(self,z,weights,args{:});
end
