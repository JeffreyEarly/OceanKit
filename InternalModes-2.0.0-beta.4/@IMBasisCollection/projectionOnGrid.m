function projection = projectionOnGrid(self,z,weights,options)
% Construct a projection on the supplied fixed grid and quadrature.
%
% Preserve explicit columns, physical signed pairing, and reference provenance.
% No quadrature fitting, acceptance policy, or mode selection is performed.
% - Topic: Construct projections
% - Parameter z: physical sample column
% - Parameter weights: quadrature weights aligned with z
% - Parameter options.page: requested page index
% - Parameter options.variable: sampled variable, formulation by default
% - Parameter options.columns: explicit array columns, all by default
% - Returns projection: fixed numerical projection
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string = ""
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
end
arguments (Output)
    projection (1,1) IMProjection
end
args = namedargs2cell(options);
projection = self.prepareProjection(z,weights,args{:});
end
