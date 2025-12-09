% Demonstrates how to use the PoissonFlowFromFlux in the WVDiagnostics

basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=9; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

energy_fluxes = wvd.exactEnergyFluxesTemporalAverage(timeIndices=51:251);
enstrophy_fluxes = wvd.exactEnstrophyFluxesTemporalAverage(timeIndices=51:251);

inertial_fluxes = wvd.quadraticEnergyTriadFluxesTemporalAverage(timeIndices=51:251);
ggg = inertial_fluxes(1).te_gmda;
wgg = inertial_fluxes(1).te_wave + inertial_fluxes(2).te_gmda + inertial_fluxes(3).te_gmda;
wwg = inertial_fluxes(2).te_wave + inertial_fluxes(3).te_wave + inertial_fluxes(4).te_gmda;
www = inertial_fluxes(4).te_wave;

%% Choose which flux to plot, and start by using a linear axes
flux = energy_fluxes(1).te/wvd.flux_scale;
flux = enstrophy_fluxes(1).Z0/wvd.z_flux_scale;
flux = www/wvd.flux_scale;

[X,Y,U,V] = wvd.PoissonFlowFromFlux(flux.');

kRadial = wvd.kRadial;
jWavenumber = wvd.jWavenumber;
figure
tiledlayout(2,2)
nexttile

pcolor(kRadial,jWavenumber,flux); shading flat;
xlim([kRadial(1) 1.6e-3])
ylim([jWavenumber(1) 1.6e-3])

colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(flux(:)))*[-1 1]/2)
colorbar("eastoutside")
hold on,
quiver(X,Y,U,V,Color=0*[1 1 1])

% Now use log axes

[logX,logY,Uprime,Vprime] = wvd.RescalePoissonFlowFluxForLogSpace(X,Y,U,V);

% figure
% tiledlayout(1,2)
nexttile
scale = 2;
quiver(logX,logY,scale*Uprime,scale*Vprime,"off",Color=0*[1 1 1]);
nexttile
quiver(X,Y,U,V,Color=0*[1 1 1])

% Now for something completely different
N = 30;
xLinLog = linspace(min(logX(:)),max(logX(:)),N);
yLinLog = linspace(min(logY(:)),max(logY(:)),N/2);
[XLinLog,YLinLog] = ndgrid(xLinLog,yLinLog);
fluxLinLog = interpn(logX,logY,(flux.'),XLinLog,YLinLog);
ULinLog = interpn(logX,logY,Uprime,XLinLog,YLinLog);
VLinLog = interpn(logX,logY,Vprime,XLinLog,YLinLog);

% figure
nexttile
pcolor(XLinLog,YLinLog,fluxLinLog); shading flat; colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(fluxLinLog(:)))*[-1 1]/2)
hold on
quiver(XLinLog,YLinLog,ULinLog,VLinLog,Color=0*[1 1 1])


%%
figure
plot(wvd.kRadial,cumsum(squeeze(sum(flux,1)))), hold on
plot(wvd.kRadial,-(squeeze(sum(U/wvd.kRadial(2),2))))

%%
figure
plot(wvd.jWavenumber,cumsum(squeeze(sum(flux,2)))), hold on
plot(wvd.jWavenumber,-(squeeze(sum(V/wvd.jWavenumber(2),1))))


%%

[X,Y,U,V] = wvd.PoissonFlowFromFluxType1(flux.');
figure
quiver(X,Y,U,V,Color=0*[1 1 1])