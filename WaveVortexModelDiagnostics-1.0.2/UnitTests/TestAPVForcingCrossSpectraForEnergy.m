% basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = "/Volumes/Samsung_T7/CimRuns_June2025/output/";
% basedir = '/Users/cwortham/Documents/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns/output/';
% basedir = '/Volumes/SanDiskExtremePro/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns_June2025_v2/output/';

runNumber=18; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");
% wvd = WVDiagnostics(basedir + getRunParameters(runNumber) + ".nc");

% wvt = wvd.wvt;
% ncfile = wvd.wvfile;

options.shouldMeasureAntialiasingFlux = true;

if options.shouldMeasureAntialiasingFlux
    [wvt_lowres, ncfile] = WVTransform.waveVortexTransformFromFile(wvd.wvfile.path,iTime=Inf);
    wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
else
    [wvt, ncfile] = WVTransform.waveVortexTransformFromFile(wvd.wvfile.path,iTime=Inf);
end

wvt.addOperation(EtaTrueOperation());
wvt.addOperation(APEOperation(wvt));
wvt.addOperation(APVOperation());
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

forcingNames = wvt.forcingNames;

%%
indices = 983;
Z2_qgpv_t = zeros(length(indices),1);
Z2_t = zeros(length(indices),1);
iIndex = 1;

if options.shouldMeasureAntialiasingFlux
    wvt_lowres.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
    wvt.t = wvt_lowres.t;
    [wvt.A0,wvt.Ap,wvt.Am] = wvt_lowres.spectralVariableWithResolution(wvt,wvt_lowres.A0,wvt_lowres.Ap,wvt_lowres.Am);
else
    wvt.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
end

% wvt.initFromNetCDFFile(ncfile,iTime=indices(iIndex));
F = wvt.fluxForForcing();

figure; tl = tiledlayout(2,1,TileSpacing="compact");
fprintf("Energy forcing:\n");
energyScale = wvd.flux_scale;
eta_true = wvt.eta_true;
totalQuadratic = 0;
totalExact = 0;
for iForce=1:length(forcingNames)

    [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(F{forcingNames(iForce)}.Fp,F{forcingNames(iForce)}.Fm,F{forcingNames(iForce)}.F0);
    if isa(wvt,"WVTransformHydrostatic")
        [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        F_density = wvt.u .* Fu + wvt.v .* Fv+ wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;

        S_energy = wvt.crossSpectrumWithFgTransform(wvt.u,Fu);
        S_energy = S_energy + wvt.crossSpectrumWithFgTransform(wvt.v,Fv);
        S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.eta_true,Feta);
    else
        [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        F_density = wvt.u .* Fu + wvt.v .* Fv +  wvt.w .* Fw + wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;

        S_energy = wvt.crossSpectrumWithFgTransform(wvt.u,Fu);
        S_energy = S_energy + wvt.crossSpectrumWithFgTransform(wvt.v,Fv);
        S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.w./wvt.N2Function(wvt.Z),Fw);
        S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.eta_true,Feta);
    end
    if forcingNames(iForce) == "nonlinear advection"
        F_density = F_density + wvt.w .* shiftdim(wvt.N2,-2) .* (wvt.eta_true-wvt.eta);
        S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.w,wvt.eta_true-wvt.eta);
    end

    E_jk = wvt.transformToRadialWavenumber(Ep+Em+E0);
    if isscalar(indices)
        totalQuadratic = totalQuadratic + sum(E_jk(:))/energyScale;
        fprintf("Total quadratic energy for " + forcingNames(iForce) + " forcing: " + sum(E_jk(:))/energyScale + "\n");
    end

    E = int_vol(F_density);
    if isscalar(indices)
        totalExact = totalExact + E/energyScale;
        fprintf("Total nonlinear energy for " + forcingNames(iForce) + " forcing: " + E/energyScale + " (" + sum(S_energy(:))/energyScale + ")\n");
    end

    nexttile(tl,1)
    E_k = wvt.transformToPseudoRadialWavenumberA0(E_jk);
    if iForce == 1
        plot(wvt.kPseudoRadial,cumsum(E_k)/energyScale,Color=0*[1 1 1],LineWidth=2), hold on
    else
        plot(wvt.kPseudoRadial,cumsum(E_k)/energyScale)
    end

    nexttile(tl,2)
    % S_f_R = wvt.transformToRadialWavenumber(S_f);
    S_jk = wvt.transformToRadialWavenumber( S_energy);
    S_k = wvt.transformToPseudoRadialWavenumberA0( S_jk );
    if iForce == 1
        plot(wvt.kPseudoRadial,cumsum(S_k)/energyScale,Color=0*[1 1 1],LineWidth=2), hold on
        dE = E_k-S_k;
    else
        plot(wvt.kPseudoRadial,cumsum(S_k)/energyScale)
    end
    

    % 
    % Z2_qgpv_t(iIndex) = sum(Z(:))/energyScale;
    % 
    % nexttile(tl,1)
    % % S_f_R = wvt.transformToRadialWavenumber(Z);
    

    % F_pv = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(Feta);
    % Z2 = wvt.qgpv .* F_pv;
    % if isscalar(indices)
    %     fprintf("Total approximate enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
    % end
    % F0_pv = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(F_pv));


    % zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
    % zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
    % zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y
    % 
    % if forcingNames(iForce) == "nonlinear advection"
    %     G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*(Feta + wvt.w);
    %     % G_eta = (-wvt.u .* wvt.diffX(eta_true) - wvt.v .* wvt.diffY(eta_true) - wvt.w .* (wvt.diffZG(eta_true) - 1));
    % else
    %     G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*Feta;
    % end
    % % G_eta = Feta;
    % FZ_L = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(G_eta);
    % FZ_NL = - wvt.zeta_x .* wvt.diffX(G_eta) - wvt.zeta_y .* wvt.diffY(G_eta)- wvt.zeta_z .* wvt.diffZG(G_eta);
    % FZ_NL = FZ_NL - wvt.diffX(eta_true) .* DF_x - wvt.diffY(eta_true) .* DF_y - wvt.diffZG(eta_true) .* DF_z;
    % FZ = FZ_L + FZ_NL;
    % Z2 = wvt.apv .* FZ;
    % 
    % Z2_t(iIndex) = int_vol(Z2)/energyScale;
    % 
    % % S_f = wvd.crossSpectrumWithFgTransform(FZ,wvt.apv);
    % 
    % phi_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.apv));
    % gamma_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(FZ));
    % S_f = prefactor .* real(phi_bar .* conj(gamma_bar));
    % 
    % nexttile(tl,2)
    % % S_f_R = wvt.transformToRadialWavenumber(S_f);
    % S_f_jK = wvt.transformToRadialWavenumber( S_f);
    % S_f_R = wvt.transformToPseudoRadialWavenumberA0( S_f_jK );
    % if iForce == 1
    %     plot(wvt.kPseudoRadial,cumsum(S_f_R)/energyScale,Color=0*[1 1 1],LineWidth=2), hold on
    % else
    %     plot(wvt.kPseudoRadial,cumsum(S_f_R)/energyScale)
    % end
    % 
    % if isscalar(indices)
    %     fprintf("Total nonlinear enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/energyScale + " (" + sum(S_f(:))/energyScale + ")\n");
    % end
    % 
end

fprintf("Total quadratic energy for all forcing: " + totalQuadratic + "\n");
fprintf("Total nonlinear energy for all forcing: " + totalExact + "\n");

nexttile(tl,1);
plot(wvt.kPseudoRadial,zeros(size(wvt.kPseudoRadial)),Color=0*[1 1 1],LineWidth=1)
xlog
nexttile(tl,2);
plot(wvt.kPseudoRadial,zeros(size(wvt.kPseudoRadial)),Color=0*[1 1 1],LineWidth=1)
plot(wvt.kPseudoRadial,dE/energyScale,Color=0*[1 1 1],LineWidth=1,LineStyle="--")
xlog
legend(forcingNames)

return
%%

flux = S_jk;

jWavenumber = 1./sqrt(wvt.Lr2);
jWavenumber(1) = 0; % barotropic mode is a mean?
[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');

figure, jpcolor(wvt.kRadial,jWavenumber,flux); shading flat;
colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(S_jk(:)))*[-1 1])
colorbar("eastoutside")
hold on,
quiver(X,Y,10*U,10*V,Color=0*[1 1 1])
xlim([0 5e-4])
ylim([0 5e-4])

%%
figure
flux = wvt.transformToRadialWavenumber(Z);
jpcolor(wvt.kRadial,jWavenumber,flux); shading flat;
colormap(WVDiagnostics.crameri('-bam'))
clim(max(abs(S_jk(:)))*[-1 1])
colorbar("eastoutside")
xlim([0 5e-4])
ylim([0 5e-4])

%%

radialWavelength = 2*pi./wvt.kRadial/1000;
radialWavelength(1) = 1.5*radialWavelength(2);

verticalWavelength = 2*pi./jWavenumber/1000;
verticalWavelength(1) = 1.5*verticalWavelength(2);

figure
ax1 =axes; 
pcolor(ax1,radialWavelength,verticalWavelength,flux); shading flat;
xscale("log"), yscale("log")
set(gca,'XDir','reverse'), set(gca,'YDir','reverse')
colormap(WVDiagnostics.crameri('-bam')), 
clim(max(abs(S_jk(:)))*[-1 1])



%%

% check K only transform
a = squeeze(mean(mean(wvt.apv .* wvt.apv,1),2));
b = squeeze(sum(sum(fft2(wvt.apv) .* conj(fft2(wvt.apv)),1),2))/(wvt.Nx*wvt.Ny)^2;
c = sum(prefactorK.*(wvt.transformFromSpatialDomainWithFourier(wvt.apv).*conj(wvt.transformFromSpatialDomainWithFourier(wvt.apv))),2);

sum(wvt.z_int .* b)

% check J only transform
apv_zxy = reshape(shiftdim(wvt.apv,2),wvt.Nz,[]);
aa = prefactorJ.* mean( ((wvt.PF0*apv_zxy)./wvt.P0).^2 ,2);
sum(aa)

% check combined transform
phi_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.apv));
S_f = prefactor .* real(phi_bar .* conj(phi_bar));
sum(S_f(:))

S_f = wvd.crossSpectrumWithFgTransform(wvt.apv,wvt.apv);
sum(S_f(:))