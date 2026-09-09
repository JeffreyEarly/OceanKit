%% Compose APV and boundary modes at positive horizontal wavenumber
% A geostrophic state at positive horizontal wavenumber has two vertical
% pieces. APV modes describe the interior field, while one canonical
% zero-APV coordinate carries each active endpoint anomaly. This example
% builds both pieces on a shared exponential-stratification grid and
% combines them with IMGeostrophicTransform.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Build the generalized-energy APV transform
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
f0 = 1e-4;
g0 = 0.03;
gd = 0.01;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nModes = 12;
nEVP = 192;
solver = IMSolverSpectral(nEVP=nEVP);
evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd,surfaceBoundary="freeSurface");
basisSet = solver.solveEVP(evp,nModes=nModes+3);

z = linspace(zDomain(1),zDomain(2),97).';
apvTransform = basisSet.discreteTransform(z=z,nModes=nModes,variables=["F","G"]);

%% Solve the canonical boundary modes on each wavenumber page
% Finite g0 and gd activate both endpoints. The zero-APV basis therefore
% contains a surface-response column and a bottom-response column for each
% positive horizontal wavenumber.
k = [8e-6 1.6e-5 3.2e-5];
zeroProblem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=zDomain,f0=f0,g=g,k=k,endpoints=["surface","bottom"],surfaceBoundary="freeSurface");
zeroAPVModes = solver.solveGeostrophicZeroAPVModes(zeroProblem);

transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd);

fprintf("\nThe transform has %d APV modes, %d wavenumber pages, and active endpoints %s.\n", ...
    nModes,numel(transform.k),join(transform.activeEndpoints,", "));
fprintf("Minimum relative mu separation: %.3e.\n",transform.compatibilityDiagnostics.minimumRelativeMuSeparation);

%% Transform an admissible state
% The second array dimension is the wavenumber page. The third dimension
% below represents two independent fields and is preserved by every
% operation.
nFields = 2;
APVCoefficients = zeros(nModes,numel(k),nFields);
for iField = 1:nFields
    for iK = 1:numel(k)
        APVCoefficients(:,iK,iField) = exp(1i*(0.2*iField+0.15*iK))*(0.72.^(0:nModes-1)).';
    end
end
zeroAPVCoefficients = reshape([0.7 -0.3 0.45 -0.18 0.25 -0.1 0.5 -0.2 0.3 -0.12 0.15 -0.06],2,numel(k),nFields);

[APV,endpointAnomalies] = transform.transformStateBack(APVCoefficients=APVCoefficients,zeroAPVCoefficients=zeroAPVCoefficients);
[APVCoefficientsBack,zeroAPVCoefficientsBack] = transform.transformStateForward(APV=APV,endpointAnomalies=endpointAnomalies);

APVError = norm(APVCoefficientsBack(:)-APVCoefficients(:))/norm(APVCoefficients(:));
zeroAPVError = norm(zeroAPVCoefficientsBack(:)-zeroAPVCoefficients(:))/norm(zeroAPVCoefficients(:));
fprintf("APV coefficient round-trip error: %.3e.\n",APVError);
fprintf("Zero-APV coefficient round-trip error: %.3e.\n",zeroAPVError);

figure(Name="Geostrophic state composition",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(real(APV(:,1,1)),z,LineWidth=1.5)
grid on
xlabel("Re(q)")
ylabel("z (m)")
title("Interior APV, first page")

nexttile
imagesc(1:nModes,k,abs(squeeze(transform.apvEndpointResponse(1,:,:))).')
axis xy
colorbar
xlabel("APV mode column")
ylabel("\kappa (m^{-1})")
title("Surface response magnitude")

%% Use rotated zero-APV coordinates
% Boundary-normalized coordinates make endpoint inversion direct. A
% boundary-depth rotation gives a different coefficient representation of
% exactly the same physical state.
depthModes = zeroAPVModes.rotateBoundaryDepth(g0=g0,gd=gd);
[~,depthCoefficients] = transform.transformStateForward(APV=APV,endpointAnomalies=endpointAnomalies,zeroAPVCoordinates=depthModes);
[APVFromDepth,endpointsFromDepth] = transform.transformStateBack(APVCoefficients=APVCoefficientsBack,zeroAPVCoefficients=depthCoefficients,zeroAPVCoordinates=depthModes);

fprintf("Rotated-coordinate APV reconstruction error: %.3e.\n",norm(APVFromDepth(:)-APV(:))/norm(APV(:)));
fprintf("Rotated-coordinate endpoint reconstruction error: %.3e.\n",norm(endpointsFromDepth(:)-endpointAnomalies(:))/norm(endpointAnomalies(:)));

%% Project a generic source
% A generic source need not satisfy the state endpoint relation. The APV
% source uses volume-plus-endpoint modal pairings, and the zero-APV source
% uses a small generalized-energy solve on each wavenumber page.
x = (z-zDomain(1))/D;
vorticitySource = zeros(length(z),numel(k),nFields);
displacementSource = zeros(length(z),numel(k),nFields);
for iField = 1:nFields
    for iK = 1:numel(k)
        phase = exp(1i*(0.25*iField+0.1*iK));
        vorticitySource(:,iK,iField) = phase*(sin(pi*iField*x)+0.2*cos(2*pi*x));
        displacementSource(:,iK,iField) = conj(phase)*(0.4+cos(pi*x)+0.1*iK*sin(3*pi*x));
    end
end

[APVSourceCoefficients,zeroAPVSourceCoefficients] = transform.transformSourceForward(vorticitySource=vorticitySource,displacementSource=displacementSource);
[APVSourceInDepthCoordinates,zeroAPVSourceInDepthCoordinates] = transform.transformSourceForward( ...
    vorticitySource=vorticitySource,displacementSource=displacementSource,zeroAPVCoordinates=depthModes);

fprintf("Generic source arrays have sizes %s and %s.\n",mat2str(size(APVSourceCoefficients)),mat2str(size(zeroAPVSourceCoefficients)));
fprintf("The APV source is coordinate invariant: %.3e.\n",norm(APVSourceInDepthCoordinates(:)-APVSourceCoefficients(:)));

figure(Name="Geostrophic source coefficients",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
semilogy(0:nModes-1,abs(APVSourceCoefficients(:,1,1)),"o-",LineWidth=1.3)
grid on
xlabel("APV mode number")
ylabel("|S_q^j|")
title("Interior source coefficients")

nexttile
bar(categorical(transform.activeEndpoints),abs(zeroAPVSourceCoefficients(:,1,1)))
grid on
ylabel("|S_0|")
title("Boundary-normalized source")
