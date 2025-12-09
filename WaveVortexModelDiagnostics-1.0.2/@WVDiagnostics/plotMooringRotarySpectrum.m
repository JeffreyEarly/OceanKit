function fig = plotMooringRotarySpectrum(self,options)
% Plot rotary spectra from mooring velocity time series.
%
% Plot rotary spectra from mooring velocity time series.
% This function assuminges mooring horizontal velocity time series are saved to file.
% This function uses mspec from jLab.
% Compute and plot negative and positive rotary spectra from mooring
% horizontal velocity records stored in the WVModel output. Reads mooring
% time and velocity variables from self.wvfile, computes multitaper rotary
% spectra for each depth, averages over mooring IDs, and draws a two-panel
% log-log figure of negative and positive rotary spectra. Annotates plots
% with inertial frequency (f), M2 tidal frequency, f+M2, and min/max buoyancy
% frequency (N). Legend, panel titles, and overall title are optional.
%
% - Topic: Figures — Ancillary
% - Declaration: fig = plotMooringRotarySpectrum(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter title: (optional) Overall figure title (default: "total velocity"; set to "none" to suppress)
% - Parameter shouldShowLegend: (optional) (optional, logical) Show legend (default: true)
% - Parameter shouldShowSpectralTitles: (optional) (optional, logical) Show per-panel titles (default: true)
% - Returns fig: Handle to the generated figure
arguments
    self WVDiagnostics
    options.visible = "on"
    options.title = "total velocity";
    options.shouldShowLegend = true;
    options.shouldShowSpectralTitles = true
end

figTitle = "total velocity";


t = self.wvfile.readVariables('mooring/t');
z = self.wvfile.readVariables('mooring/mooring_z');
u = self.wvfile.readVariables('mooring/mooring_u');
v = self.wvfile.readVariables('mooring/mooring_v');

cmap = self.cmocean('deep',length(z));
wvt = self.wvt;
legendCell = cellstr(strcat(string(round(abs(flipud(z)))),' m'));
M2Period = 12.420602*3600; % M2 tidal period, s
purple = [0.4940, 0.1840, 0.5560]; % extra color for labels

% compute rotary spectrum at each (mooring_z,mooring_id)
% cv_mooring = shiftdim( u + sqrt(-1)*v, 2);
cv_mooring = shiftdim( u + sqrt(-1)*v, 2);
[PSI,LAMBDA] = sleptap(length(t));
[omega_p, Spp, Snn, Spn] = mspec(t(2)-t(1),cv_mooring,PSI);
% omega = [ -flipud(omega_p(2:end)); omega_p];
% average over mooring_id
Spp = mean(Spp,3);
Snn = mean(Snn,3);
% Spn = mean(Spn,3);
% We want the integral of this to give us the variance back, so we need to
% divide by 2*pi
% Srot = (1/(2*pi))*[flipud(Snn); Spp(2:end,:)];

% two panel rotary spectrum
fig = figure(Visible=options.visible);
tl = tiledlayout(1,2,"TileSpacing","tight");
if options.title ~= "none"
    title(tl, options.title, 'Interpreter', 'none')
end
ax1 = nexttile(tl,1);
set(ax1,'xdir','reverse')
set(ax1,'XScale','log')
set(ax1,'YScale','log')
axis tight
box on
hold on
if options.shouldShowSpectralTitles
title('Negative Rotary Spectra')
end
ax2 = nexttile(tl,2);
set(ax2,'XScale','log')
set(ax2,'YScale','log')
axis tight
box on
hold on
if options.shouldShowSpectralTitles
title('Positive Rotary Spectra')
end

% plot model spectrum
set(ax1,'ColorOrderIndex',1)
set(ax2,'ColorOrderIndex',1)
set(ax1,'ColorOrder',cmap)
set(ax2,'ColorOrder',cmap)
% negative
plt1 = loglog(ax1,omega_p*86400/2/pi,fliplr(Snn),'linewidth',1);
set(plt1,{'DisplayName'},legendCell)
% positive
plt2 = loglog(ax2,omega_p*86400/2/pi,fliplr(Spp),'linewidth',1);
set(plt2,{'DisplayName'},legendCell)
% add f, M2, N lines
plt4 = plot(ax1,[wvt.f,wvt.f]*86400/2/pi,ylim,'r:','DisplayName','f');
plt5 = plot(ax2,[wvt.f,wvt.f]*86400/2/pi,ylim,'r:','DisplayName','f');
plt6 = plot(ax1,[2*pi/M2Period,2*pi/M2Period]*86400/2/pi,ylim,'b:','DisplayName','M2');
plt7 = plot(ax2,[2*pi/M2Period,2*pi/M2Period]*86400/2/pi,ylim,'b:','DisplayName','M2');
plt8 = plot(ax1,[2*pi/M2Period+wvt.f,2*pi/M2Period+wvt.f]*86400/2/pi,ylim,'Color',purple,'linestyle',':','DisplayName','f+M2');
plt9 = plot(ax2,[2*pi/M2Period+wvt.f,2*pi/M2Period+wvt.f]*86400/2/pi,ylim,'Color',purple,'linestyle',':','DisplayName','f+M2');
plt10 = plot(ax1,[min(sqrt(wvt.N2)),min(sqrt(wvt.N2))]*86400/2/pi,ylim,'k:','DisplayName','min N(z)');
plt11 = plot(ax2,[min(sqrt(wvt.N2)),min(sqrt(wvt.N2))]*86400/2/pi,ylim,'k:','DisplayName','min N(z)');
plt12 = plot(ax1,[max(sqrt(wvt.N2)),max(sqrt(wvt.N2))]*86400/2/pi,ylim,'k:','DisplayName','max N(z)');
plt13 = plot(ax2,[max(sqrt(wvt.N2)),max(sqrt(wvt.N2))]*86400/2/pi,ylim,'k:','DisplayName','max N(z)');
% power law reference
plt14 = loglog(ax1,omega_p*86400/2/pi, 10*(omega_p*86400/2/pi).^-2 ,'k:','LineWidth',1,'DisplayName','\omega^{-2}');
plt15 = loglog(ax2,omega_p*86400/2/pi, 10*(omega_p*86400/2/pi).^-2 ,'k:','LineWidth',1,'DisplayName','\omega^{-2}');
% tweak plot
xlabel(tl,'frequency (cycles per day)')
ylabel(tl,'Horizontal velocity spectrum (m^2/s^2)')
ax1.YLim = ([min([ax1.YLim,ax2.YLim]),max([ax1.YLim,ax2.YLim])]);
ax2.YLim = ax1.YLim;
ax1.XLim = ([omega_p(2)*86400/2/pi,1.2*max([ax1.XLim,ax2.XLim])]);
ax2.XLim = ax1.XLim;
ax2.YTickLabels={};
% legend and labels
% legend(ax2,[plt2',plt17,plt5,plt7,plt9,plt11,plt13,plt15],'location','northeastoutside')
if options.shouldShowLegend
    legend(ax2,plt2','location','eastoutside')
end
textY = .5*max(ax2.YLim);
text(ax2,.9*wvt.f*86400/2/pi,textY,'f','Color','r')
text(ax2,.75*2*pi/M2Period*86400/2/pi,textY,'M2','Color','b')
text(ax2,.7*(2*pi/M2Period+wvt.f)*86400/2/pi,textY,'f+M2','Color',purple)
text(ax2,1.05*min(sqrt(wvt.N2))*86400/2/pi,textY,'min N(z)','Color','k')
text(ax2,.5*max(sqrt(wvt.N2))*86400/2/pi,textY,'max N(z)','Color','k')
text(ax2,1.1*omega_p(3)*86400/2/pi, 10*(omega_p(3)*86400/2/pi).^-2,'\omega^{-2}','Color','k')
text(ax1,.9*wvt.f*86400/2/pi,textY,'f','Color','r','HorizontalAlignment','right')
text(ax1,.75*2*pi/M2Period*86400/2/pi,textY,'M2','Color','b','HorizontalAlignment','right')
text(ax1,.7*(2*pi/M2Period+wvt.f)*86400/2/pi,textY,'f+M2','Color',purple,'HorizontalAlignment','right')
text(ax1,1.05*min(sqrt(wvt.N2))*86400/2/pi,textY,'min N(z)','Color','k','HorizontalAlignment','right')
text(ax1,.5*max(sqrt(wvt.N2))*86400/2/pi,textY,'max N(z)','Color','k','HorizontalAlignment','right')
text(ax1,1.1*omega_p(3)*86400/2/pi, 10*(omega_p(3)*86400/2/pi).^-2,'\omega^{-2}','Color','k','HorizontalAlignment','right')

end