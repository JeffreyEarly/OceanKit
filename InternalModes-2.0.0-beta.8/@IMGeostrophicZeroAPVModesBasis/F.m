function values = F(self,z,options)
% Evaluate the aligned F boundary-buoyancy profile.
% - Topic: Evaluate boundary-buoyancy modes
% - Parameter z: physical sample column
% - Parameter options.pages: source wavenumber pages in caller order
% - Parameter options.normalization: optional per-call normalization
% - Returns values: sample-by-column-by-page values
arguments
    self IMGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
    options.normalization (1,1) string = self.normalization
end
values = self.evaluate(z,variable="F",derivativeOrder=0,pages=options.pages,normalization=options.normalization);
end
