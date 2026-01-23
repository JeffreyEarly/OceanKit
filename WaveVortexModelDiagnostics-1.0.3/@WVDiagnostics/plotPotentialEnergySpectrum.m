function fig = plotPotentialEnergySpectrum(self,options)
% Plot the potential energy spectrum
%
% Attempts to plot the potential energy spectrum at a given time index, although this might not be fully implemented yet.
%
% - Topic: Figures — Model Snapshot
% - Declaration: fig = plotEnergyOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter iTime: time-related parameter `iTime`
% - Parameter visible: (optional) figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.iTime
    options.visible = "on"
end

if isfield(options,"iTime")
    self.iTime = options.iTime;
end

wvt = self.wvt;

% PE
PE_j_kl = 0.5*wvt.spectrumWithGgTransform(wvt.eta);
PE_j_kl(:,1) = 0;
PE_j_kR = wvt.transformToRadialWavenumber(PE_j_kl);
PE_kR = sum(PE_j_kR,1);
PE_j = sum(PE_j_kR,2);

% APE
rho_nm = wvt.chebfunForZArray(wvt.rho_nm)/wvt.rho0 - 1;
p_nm = - wvt.g * cumsum(rho_nm);
p_nm = p_nm - p_nm(0);
first_term = wvt.crossSpectrumWithGgTransform(wvt.g*wvt.eta_true./wvt.N2Function(wvt.Z),rho_nm(wvt.Z - wvt.eta_true));
second_term = wvt.crossSpectrumWithFgTransform(p_nm(wvt.Z) - p_nm(wvt.Z - wvt.eta_true),ones(size(wvt.Z)));
APE_j_kl = first_term + second_term;
APE_j_kR = wvt.transformToRadialWavenumber(APE_j_kl);
APE_kR = sum(APE_j_kR,1);
APE_j = sum(APE_j_kR,2);

options.axes = "k-radial";
switch options.axes
    case "k-pseudo-isotropic"
        radialWavelength = 2*pi./wvt.kPseudoRadial/1000;
    otherwise
        radialWavelength = 2*pi./wvt.kRadial/1000;
end
radialWavelength(1) = 1.5*radialWavelength(2);

fig = figure(Visible=options.visible);

% plot horizontal wavenumber spectrum
% axK = nexttile(n); n = n+1;
%%
nexttile
plot(radialWavelength,APE_kR), hold on
plot(radialWavelength,PE_kR),
set(gca,'XDir','reverse')
xscale('log'); yscale('log')
axis tight
title('Radial Wavenumber Spectrum')
ylabel({'potential energy (m^3 s^{-2})'});
xlabel('wavelength (km)')
% set(gca,'XTick',[])
% ylim(10.^options.clim);
% legEntry{1} = "exact APE (" + sum(APE_j_kR(:)) + ")";
% legEntry{2} = "quadratic APE (" + sum(sum(PE_j_kR(:,2:end))) + ")";
legEntry{1} = "exact APE (" + sum(APE_kR(:)) + ")";
legEntry{2} = "quadratic APE (" + sum(PE_kR(:)) + ")";
legend(legEntry);

end