function fig = plotPoissonFlowOverPcolor(wvd,options)
arguments
    wvd WVDiagnostics
    options.visible = "on"
    options.vectorFlux
    options.flux
    options.overSaturationFactor = 2;
end

color_axis_limits = max(abs(options.flux(:)))*[-1 1]/options.overSaturationFactor;

% We will pretend the "0" wavenumber is actually evenly spaced
% from the nearest two wavenumbers
kPseudoLocation = wvd.kRadial;
kPseudoLocation(1) = exp(-log(kPseudoLocation(3)) + 2*log(kPseudoLocation(2)));
jPseudoLocation = wvd.jWavenumber;
jPseudoLocation(1) = exp(-log(jPseudoLocation(3)) + 2*log(jPseudoLocation(2)));

% For interpolation to work correctly we need to repeat the
% first entry, but properly back at zero
kPadded = cat(1,0,kPseudoLocation);
jPadded= cat(1,0,jPseudoLocation);
[KPadded,JPadded] = ndgrid(kPadded,jPadded);
fluxPadded = cat(1,options.flux(1,:),options.flux);
fluxPadded = cat(2,fluxPadded(:,1),fluxPadded);

% We will use this axis to display. Widen the box so that it is
% the same size as its neighbors.
kMin = exp(-1.5*log(wvd.kRadial(3)) + 2.5*log(wvd.kRadial(2)));
jMin = exp(-1.5*log(wvd.jWavenumber(3)) + 2.5*log(wvd.jWavenumber(2)));

N = 500;
kLinLog = linspace(log10(kMin),max(log10(wvd.kRadial(:))),N);
jLinLog = linspace(log10(jMin),max(log10(wvd.jWavenumber(:))),N/2);
[KLinLog,JLinLog] = ndgrid(kLinLog,jLinLog);
fluxLinLog = interpn(KPadded,JPadded,(fluxPadded.'),10.^KLinLog,10.^JLinLog,"nearest");

fig = figure(Units='points',Visible = options.visible);
set(gcf,'PaperPositionMode','auto')

pcolor(KLinLog,JLinLog,fluxLinLog); shading flat;
colormap(WVDiagnostics.crameri('-bam'))
clim(color_axis_limits)
colorbar("eastoutside")

[X,Y,U,V] = wvd.PoissonFlowFromFlux(options.vectorFlux.');
[logX,logY,Uprime,Vprime] = wvd.RescalePoissonFlowFluxForLogSpace(X,Y,U,V,shouldOnlyRescaleDirection=false);
% we need two adjustments. First, we need to move the first row and column
% half an increment
logX(1,:) = log10(kPseudoLocation(1));
logY(:,1) = log10(jPseudoLocation(1));
scale = 1.5;
hold on
% quiver(logX,logY,scale*Uprime,scale*Vprime,'off',Color=0*[1 1 1],AutoScale="off",LineWidth=2)
quiver(logX,logY,Uprime,Vprime,Color=0*[1 1 1],AutoScale="off",LineWidth=2)
end