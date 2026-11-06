function basis = createIMTestBasis(family,background,options)
% Construct a test fixture through the public capable-family API.
arguments
    family (1,1) string
    background (1,1) IMBackground
    options.solver = IMSolverSpectral(nEVP=64)
    options.nModes (1,1) double = 3
    options.k (1,1) double = 1e-3
    options.omega (1,1) double = .005
    options.energyBoundary = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf)
    options.normalization (1,1) string = ""
    options.intervals (:,2) double = [-1000 0]
end
switch family
    case "wave"
        problem = IMInternalModes.waveModesAtWavenumber(background,k=options.k);
    case "frequency"
        problem = IMInternalModes.waveModesAtFrequency(background,omega=options.omega);
    case "inertial"
        problem = IMInternalModes.waveModesAtWavenumber(background,k=0);
    case "hydrostaticF"
        problem = IMInternalModes.hydrostaticFModes(background);
    case "hydrostaticG"
        problem = IMInternalModes.hydrostaticGModes(background);
    case "APV"
        problem = IMInternalModes.geostrophicAPVModes(background,g0=options.energyBoundary.g0,gd=options.energyBoundary.gd);
    case "MDA"
        problem = IMInternalModes.meanDensityAnomalyModes(background,g0=options.energyBoundary.g0,gd=options.energyBoundary.gd);
    case "enstrophy"
        problem = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(background,k=options.k);
    case "neutral"
        basis = IMNeutralLabelModesBasis(background,neutralIntervals=options.intervals,nModesPerInterval=options.nModes,g0=options.energyBoundary.g0,gd=options.energyBoundary.gd);
    otherwise
        error("InternalModesTests:UnknownFixture","Unknown family %s.",family);
end
if family ~= "neutral", basis = options.solver.solveEVP(problem,nModes=options.nModes); end
if options.normalization ~= "", basis = basis.withNormalization(options.normalization); end
end
