function S_f = spectrumWithFgTransform(self,f,options)
% Spectrum With Fg Transform.
%
% spectrumWithFgTransform is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Computing Spectra
% - Declaration: S_f = spectrumWithFgTransform(self,f,options)
% - Parameter self: WVDiagnostics object
% - Parameter f: input argument `f`
% - Parameter useExplicitAntialiasedWVT: (optional) input argument `useExplicitAntialiasedWVT` (default: false)
% - Returns S_f: output value `S_f`
arguments
    self WVDiagnostics
    f
    options.useExplicitAntialiasedWVT = false
end
if options.useExplicitAntialiasedWVT
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end
prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(f));
S_f = prefactor.*abs(f_bar).^2;
end