basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=1; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
energy_fluxes = wvd.exactEnergyFluxesTemporalAverage(timeIndices=51:251);


%%
flux = energy_fluxes(1).te/wvd.flux_scale;

jWavenumber = 1./sqrt(wvd.Lr2);
jWavelength = 2*pi./wvd.jWavenumber/1000;
jWavelength(1) = 1.5*jWavelength(2);

radialWavelength = 2*pi./wvd.kRadial/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFluxDCTI(wvd.kRadial,wvd.jWavenumber,flux.');

kRadial = wvd.kRadial;
kRadial = kRadial + (kRadial(2)-kRadial(1))/2;
% kRadial(1) = 0.5*wvd.kRadial(2);
jWavenumber = wvd.jWavenumber;
jWavenumber = jWavenumber + (jWavenumber(2)-jWavenumber(1))/2;

% jWavenumber(1) = 0.5*wvd.jWavenumber(2);

figure, pcolor(kRadial,jWavenumber,flux); shading flat;

% pbaspect([wvd.jWavenumber(2)/wvd.kRadial(2) 1 1])

% set(gca,'YDir','reverse')
% set(gca,'XDir','reverse')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
xlim([kRadial(1) 1.6e-3])
ylim([jWavenumber(1) 1.6e-3])

colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(flux(:)))*[-1 1]/100)
colorbar("eastoutside")
hold on,
quiver(X,Y,U,V,Color=0*[1 1 1])

%%
Xaxes = X;
Xaxes(1,:) = kRadial(1);
Yaxes = Y;
Yaxes(:,1) = jWavenumber(1);
Xaxes = log(Xaxes);
Yaxes = log(Yaxes);

figure
scale = 1e6;
quiver(Xaxes,Yaxes,scale*U./abs(Xaxes),scale*V./abs(Yaxes),"off",Color=0*[1 1 1])
% quiver(Xaxes,Yaxes,scale*U./abs(Xaxes),scale*V./abs(Yaxes),Color=0*[1 1 1])
xlim([Xaxes(1,1) log(1.6e-3)])
ylim([Yaxes(1,1) log(1.6e-3)])

%%

figure, plot(kRadial,cumsum(squeeze(sum(flux,1)))), xlog
figure, plot(kRadial,(squeeze(mean(U,2)))), xlog

%%
Nj = size(flux,1);
Nk = size(flux,2);
flux_zero_pad = cat(1,zeros(1,Nk),flux,zeros(1,Nk));
flux_zero_pad = cat(2,zeros(Nj+2,1),flux_zero_pad,zeros(Nj+2,1));
dk = wvd.kRadial(2)-wvd.kRadial(1);
kRadial_zero_pad = cat(1,0,wvd.kRadial+dk,wvd.kRadial(end)+2*dk);
dj = wvd.jWavenumber(2)-wvd.jWavenumber(1);
j_zero_pad = cat(1,0,wvd.jWavenumber+dj,wvd.jWavenumber(end)+2*dj);

[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFluxDCTI(kRadial_zero_pad,j_zero_pad,flux_zero_pad.');
figure
quiver(X,Y,U,V,Color=0*[1 1 1])


% figure, plot(kRadial,cumsum(squeeze(sum(flux,1)))), xlog
figure, plot(kRadial_zero_pad,(squeeze(mean(U,2)))), xlog