function values = evaluate(self,z,options)
% Evaluate aligned boundary-buoyancy profiles and their regular derivatives.
%
% A per-call normalization preserves stored defaults and coordinate meaning.
% F supports derivative orders zero and one; G supports zero through two.
% - Topic: Evaluate boundary-buoyancy modes
% - Parameter z: physical sample column
% - Parameter options.variable: required F or G channel
% - Parameter options.derivativeOrder: nonnegative derivative order, zero by default
% - Parameter options.pages: source wavenumber pages in caller order
% - Parameter options.normalization: optional per-call profile normalization
% - Returns values: nZ-by-nColumns-by-nSelectedPages values
arguments (Input)
    self IMGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
    options.normalization (1,1) string = self.normalization
end
if options.variable == ""
    error("IMGeostrophicZeroAPVModesBasis:MissingVariable","Specify variable= followed by F or G.");
end
if ~ismember(options.variable,["F","G"])
    error("IMGeostrophicZeroAPVModesBasis:UnsupportedVariable","Supported channels are F and G.");
end
if (options.variable == "F" && options.derivativeOrder > 1) || options.derivativeOrder > 2
    error("IMGeostrophicZeroAPVModesBasis:UnsupportedDerivative","F supports orders zero and one; G supports zero through two.");
end
if any(options.pages > numel(self.k))
    error("IMGeostrophicZeroAPVModesBasis:InvalidPage","pages must select existing source wavenumber pages.");
end
if any(z<self.zDomain(1) | z>self.zDomain(2))
    error("IMGeostrophicZeroAPVModesBasis:OutsideDomain","Samples must lie inside zDomain.");
end
factors = self.normalizationFactors(options.normalization,pages=options.pages);
if options.derivativeOrder == 0
    values = self.endpointUnscaledProfiles(options.variable,z,options.pages);
elseif options.variable == "F"
    values = -self.N2(z)/self.g.*self.endpointUnscaledProfiles("G",z,options.pages);
elseif options.derivativeOrder == 1
    values = reshape(-self.g*self.k(options.pages).^2/self.f0^2,1,1,[]).*self.endpointUnscaledProfiles("F",z,options.pages);
else
    values = reshape(self.k(options.pages).^2/self.f0^2,1,1,[]).*self.N2(z).*self.endpointUnscaledProfiles("G",z,options.pages);
end
values = values./reshape(factors,1,self.nColumns,numel(options.pages));
end
