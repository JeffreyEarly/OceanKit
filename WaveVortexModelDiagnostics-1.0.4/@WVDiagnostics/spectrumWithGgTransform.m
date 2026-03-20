function S_f = spectrumWithGgTransform(self,f)
% Spectrum With Gg Transform.
%
% spectrumWithGgTransform is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Computing Spectra
% - Declaration: S_f = spectrumWithGgTransform(self,f)
% - Parameter self: WVDiagnostics object
% - Parameter f: input argument `f`
% - Returns S_f: output value `S_f`
arguments
    self WVDiagnostics
    f
end
wvt = self.wvt;

prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = wvt.g * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithGg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor .* abs(f_bar).^2;
end