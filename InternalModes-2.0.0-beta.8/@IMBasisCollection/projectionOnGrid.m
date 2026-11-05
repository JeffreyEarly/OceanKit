function projection = projectionOnGrid(self,z,weights,options)
% Construct a physical scalar projection on the supplied fixed rule.
%
% Preserve explicit columns, physical signed pairing, and reference provenance.
% No quadrature fitting, acceptance policy, or mode selection is performed.
% Neutral labels accept APV samples and use their prescribed continuous
% dual, with an identity row-scaled target; supplied quadrature errors remain.
% The selected basis must provide the requested physical scalar capability.
% Auxiliary EVP pairings remain available through collection assessment.
% - Topic: Construct projections
% - Parameter z: physical sample column
% - Parameter weights: quadrature weights aligned with z
% - Parameter options.page: requested page index
% - Parameter options.variable: required independent scalar channel
% - Parameter options.columns: explicit array columns, all by default
% - Parameter options.prepareRecovery: also prepare numerical recovery for Galerkin records
% - Returns projection: fixed numerical projection
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string = ""
    options.columns = ":"
    options.prepareRecovery (1,1) logical = false
end
arguments (Output)
    projection (1,1) IMProjection
end
if options.variable == ""
    error("IMScalarProjection:MissingVariable","Specify the independent scalar observations with variable=.");
end
basis = self.selectProjectionBasis(options.page,options.columns);
physicalOptions = rmfield(options,["page","columns"]);
args = namedargs2cell(physicalOptions);
projection = IMScalarProjectionTools.projectionOnGrid(basis,z,weights,args{:});
end
