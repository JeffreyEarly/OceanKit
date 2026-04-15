function S_f = crossSpectrumWithGgTransform(self,phi,gamma)
% Cross Spectrum With Gg Transform.
%
% crossSpectrumWithGgTransform is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Computing Spectra
% - Declaration: S_f = crossSpectrumWithGgTransform(self,phi,gamma)
% - Parameter self: WVDiagnostics object
% - Parameter phi: input argument `phi`
% - Parameter gamma: input argument `gamma`
% - Returns S_f: output value `S_f`
arguments
    self WVDiagnostics
    phi
    gamma
end
wvt = self.wvt;

prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = wvt.g * prefactorK;

phi_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(phi));
gamma_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(gamma));
S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
end