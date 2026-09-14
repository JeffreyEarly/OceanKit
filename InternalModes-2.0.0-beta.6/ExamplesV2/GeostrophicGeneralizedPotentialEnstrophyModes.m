%% Generalized-potential-enstrophy modes and grids at several wavenumbers
% At fixed horizontal wavenumber, the generalized-potential-enstrophy
% modes jointly diagonalize physical energy and the positive interior-plus-
% endpoint potential-enstrophy metric. This example follows the grid-design
% workflow in `WaveVortexVerticalGridDesign.m`, but lets each wavenumber
% design its own mode-root points and fitted quadrature weights.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Specify the stratification, horizontal scales, and point count
D = 4000;
N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g = 9.81;
Nz = 64;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);
wavelength = [500e3 100e3 20e3];
k = 2*pi./wavelength;

%% Solve each wavenumber-dependent family and fit its point rule
nAvailableModes = Nz+8;
nEVP = 256;
nK = numel(k);
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSets = cell(1,nK);
zGrid = cell(1,nK);
quadratureWeights = cell(1,nK);
gridDesigns = cell(1,nK);
weightFits = cell(1,nK);
representedModeCount = zeros(nK,1);
gramError = zeros(nK,1);
minimumWeight = zeros(nK,1);
alpha0 = zeros(nK,1);
bEffective = zeros(nK,1);

for iK = 1:nK
    evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
        N2=N2,zDomain=zDomain,k=k(iK),f0=f0,g=g);
    basisSets{iK} = solver.solveEVP(evp,nModes=nAvailableModes);
    [zGrid{iK},gridDesigns{iK}] = basisSets{iK}.modeRootGrid(nPoints=Nz);
    representedModeCount(iK) = gridDesigns{iK}.representedModeCount;
    [quadratureWeights{iK},weightFits{iK}] = basisSets{iK}.quadratureWeightsForPoints( ...
        z=zGrid{iK},nModes=representedModeCount(iK),variables="F");
    gramError(iK) = weightFits{iK}.transform.relativeGramOperatorError(variable="F");
    minimumWeight(iK) = min(quadratureWeights{iK});
    alpha0(iK) = evp.parameters.alpha0;
    bEffective(iK) = evp.parameters.bEffective;

    if string(basisSets{iK}.normalization) ~= "generalizedPotentialEnstrophy"
        error("Example:UnexpectedNormalization", "The generalized-potential-enstrophy factory must install its dedicated default normalization.");
    end
end

gridSummary = table(k(:),wavelength(:)/1000,representedModeCount,alpha0,bEffective,gramError,minimumWeight, ...
    VariableNames=["k" "wavelengthKm" "representedModeCount" "alpha0" "bEffective" "FGramError" "minimumWeight"]);
disp(gridSummary)

%% Compare the independently designed point rules
% `modeRootGrid` includes both physical endpoints and the interior roots of
% the next F mode. The fitted weights then reproduce the endpoint-inclusive
% generalized-potential-enstrophy Gram matrix for the represented band.
zPlot = linspace(zDomain(1),zDomain(2),801).';
colors = colororder;
labels = compose("%g km",wavelength/1000);

figure(Name="Generalized-potential-enstrophy vertical grids",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
plot(sqrt(N2(zPlot)),zPlot,"k-",LineWidth=1.5,DisplayName="N(z)")
hold on
for iK = 1:nK
    plot(sqrt(N2(zGrid{iK})),zGrid{iK},"o",Color=colors(iK,:), ...
        MarkerSize=4,MarkerFaceColor="w",DisplayName=labels(iK));
end
hold off
grid on
xlabel("N (s^{-1})")
ylabel("z (m)")
title(sprintf("Independent %d-point grids",Nz))
legend(Location="best")

nexttile
hold on
for iK = 1:nK
    zMidpoint = 0.5*(zGrid{iK}(1:end-1)+zGrid{iK}(2:end));
    plot(diff(zGrid{iK}),zMidpoint,"o-",Color=colors(iK,:), ...
        LineWidth=1.1,MarkerSize=3,DisplayName=labels(iK));
end
hold off
grid on
xlabel("vertical spacing, \Delta z (m)")
ylabel("z (m)")
title("Mode-root grid spacing")
legend(Location="best")

nexttile
hold on
for iK = 1:nK
    plot(quadratureWeights{iK},zGrid{iK},"o-",Color=colors(iK,:), ...
        LineWidth=1.1,MarkerSize=3,DisplayName=labels(iK));
end
hold off
grid on
xlabel("fitted weight (m)")
ylabel("z (m)")
title("F-metric quadrature weights")
legend(Location="best")

%% Plot the lowest modes with each family's quadrature points
% The stored modes satisfy W_alpha/D=1. Each displayed column is divided
% by its maximum absolute value only so all four shapes remain visible.
nModeShapes = 4;
figure(Name="Generalized-potential-enstrophy modes by wavenumber",Color="w");
tiledlayout(1,nK,TileSpacing="compact",Padding="compact");
for iK = 1:nK
    FPlot = basisSets{iK}.F(zPlot);
    FPlot = FPlot(:,1:nModeShapes);
    modeScale = max(abs(FPlot),[],1);
    FPlot = FPlot./modeScale;
    FGrid = basisSets{iK}.F(zGrid{iK});
    FGrid = FGrid(:,1:nModeShapes)./modeScale;

    nexttile
    hold on
    for iMode = 1:nModeShapes
        plot(FPlot(:,iMode),zPlot,LineWidth=1.3,Color=colors(iMode,:), ...
            DisplayName="mode "+string(basisSets{iK}.modeNumber(iMode)));
        plot(FGrid(:,iMode),zGrid{iK},"o",Color=colors(iMode,:), ...
            MarkerSize=3,MarkerFaceColor="w",HandleVisibility="off");
    end
    hold off
    grid on
    xlabel("individually scaled F_j")
    ylabel("z (m)")
    title(sprintf("wavelength %g km",wavelength(iK)/1000))
    legend(Location="best")
end
