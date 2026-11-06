function values = Gzz(self,z,options)
% Evaluate the aligned Gzz boundary-buoyancy profile.
% - Topic: Evaluate boundary-buoyancy modes
% - Parameter z: physical sample column
% - Parameter options.pages: source wavenumber pages in caller order
% - Parameter options.normalization: optional per-call normalization
% - Returns values: sample-by-column-by-page values
arguments
    self IMAnalyticalGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
    options.normalization (1,1) string = self.normalization
end
values = self.evaluate(z,variable="G",derivativeOrder=2,pages=options.pages,normalization=options.normalization);
end
