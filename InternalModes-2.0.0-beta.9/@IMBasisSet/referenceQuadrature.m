function [z,weights] = referenceQuadrature(self,zBounds,refinement)
% Integrate evaluated profiles independently of the native nodal rule.
%
% Spectral profiles use an oversampled rule in the original coordinate.
% Finite-difference profiles use a rule on every reconstruction interval;
% subclasses include companion-field knots through referenceBreakpoints.
% Refinement increases the rule order without certifying convergence for
% arbitrary coefficient functions or diagnostic profile relations.
arguments (Input)
    self (1,1) IMBasisSet
    zBounds (1,2) double {mustBeReal,mustBeFinite} = self.zDomain
    refinement (1,1) double {mustBeInteger,mustBePositive} = 1
end
arguments (Output)
    z (:,1) double
    weights (:,1) double
end
if zBounds(1) >= zBounds(2)
    error("IMBasisSet:InvalidInterval","zBounds must be increasing.");
end
if zBounds(1) < self.zDomain(1) || zBounds(2) > self.zDomain(2)
    error("IMBasisSet:OutsideDomain","Reference integration bounds must lie within zDomain.");
end
if isa(self.solver,"IMSolverSpectral")
    fullBounds = [min(self.solver.xNative),max(self.solver.xNative)];
    xBounds = self.solver.xOfZ(zBounds);
    xBounds = min(max(xBounds,fullBounds(1)),fullBounds(2));
    [x,nativeWeights] = legpts(refinement*max(128,2*self.solver.nEVP),xBounds);
    z = self.solver.zOfX(x);
    coordinate = zeros(self.solver.nEVP,1);
    coordinate(1:2) = [mean(fullBounds);diff(fullBounds)/2];
    jacobian = self.solver.evaluatePhysicalDerivative(coordinate,z,1);
    weights = nativeWeights(:)./jacobian;
elseif isa(self.solver,"IMSolverFiniteDifference")
    points = self.referenceBreakpoints(zBounds);
    [nodes,rule] = legpts(8*refinement);
    widths = diff(points).';
    z = points(1:end-1).'+((nodes(:)+1)/2).*widths;
    weights = (rule(:)/2).*widths;
    z = z(:);
    weights = weights(:);
else
    [z,weights] = self.solver.nativeQuadratureRule(zBounds);
end
end
