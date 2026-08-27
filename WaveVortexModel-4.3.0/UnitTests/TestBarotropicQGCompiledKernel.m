classdef TestBarotropicQGCompiledKernel < matlab.unittest.TestCase
    properties (SetAccess = private)
        ToleranceDump (1,1) string
        ForcingDump (1,1) string
    end

    methods (TestClassSetup)
        function buildStandaloneKernel(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            scriptPath = fullfile(repositoryRoot,"tools","compiled-kernel","run_contract_tests.sh");
            [status,output] = systemWithoutMatlabRuntime(sprintf('"%s"',scriptPath));
            testCase.assertEqual(status,0,output);

            buildDirectory = fullfile(repositoryRoot,"tools","compiled-kernel","build","portable-runtime");
            configure = "cmake -S " + shellQuote(fullfile(repositoryRoot,"PortableRuntime")) + ...
                " -B " + shellQuote(buildDirectory) + ...
                " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON";
            [status,output] = systemWithoutMatlabRuntime(configure);
            testCase.assertEqual(status,0,output)
            [status,output] = systemWithoutMatlabRuntime("cmake --build " + ...
                shellQuote(buildDirectory) + " --parallel --target WVBarotropicQGToleranceDump WVBarotropicQGForcingDump");
            testCase.assertEqual(status,0,output)
            testCase.ToleranceDump = fullfile(buildDirectory,"WVBarotropicQGToleranceDump");
            testCase.assertTrue(isfile(testCase.ToleranceDump))
            testCase.ForcingDump = fullfile(buildDirectory,"WVBarotropicQGForcingDump");
            testCase.assertTrue(isfile(testCase.ForcingDump))
        end
    end

    methods (Test,TestTags="full")
        function matlabAuthoritativeParity(testCase)
            grids = [8 6; 9 7; 8 7; 9 6];
            maximumRelativeError = 0;
            maximumToleranceRelativeError = 0;
            maximumForcingRelativeError = 0;
            maximumIntegrationRelativeError = 0;
            for iGrid = 1:size(grids,1)
                for j = [0 1]
                    for shouldAntialias = [false true]
                        definition = struct( ...
                            "Nxy",grids(iGrid,:), ...
                            "Lxy",[17e3 11e3], ...
                            "h",0.8, ...
                            "j",j, ...
                            "g",9.81, ...
                            "rotationRate",7.2921e-5, ...
                            "latitude",33, ...
                            "shouldAntialias",shouldAntialias, ...
                            "planetaryRadius",6.371e6);
                        actual = barotropicFixtureDump(definition);
                        wvt = WVTransformBarotropicQG( ...
                            definition.Lxy,definition.Nxy, ...
                            h=definition.h,j=definition.j,g=definition.g, ...
                            rotationRate=definition.rotationRate, ...
                            latitude=definition.latitude, ...
                            shouldAntialias=definition.shouldAntialias, ...
                            planetaryRadius=definition.planetaryRadius);
                        diagnostic = sprintf("%dx%d, j=%d, antialias=%d",definition.Nxy,j,shouldAntialias);

                        tolerance = barotropicToleranceDump(testCase.ToleranceDump,definition,2e-7);
                        model = WVModel(wvt,shouldUseLinearDynamics=true);
                        coefficients = WVCoefficients(model,absTolerance=tolerance.absoluteToleranceScale);
                        matlabTolerance = coefficients.absErrorTolerance();
                        matlabTolerance = matlabTolerance{1};
                        testCase.verifyEqual(tolerance.absoluteTolerance(1),1.0,"C++ constrained zero mode: " + diagnostic)
                        testCase.verifyEqual(matlabTolerance(1),1.0,"MATLAB constrained zero mode: " + diagnostic)
                        maximumToleranceRelativeError = max(maximumToleranceRelativeError,verifyElementwiseRelative(testCase,tolerance.absoluteTolerance(:),matlabTolerance(:),1e-12,"adaptive tolerance vector: " + diagnostic));

                        testCase.verifyEqual(actual.contractVersion,4,diagnostic)
                        testCase.verifyEqual(actual.Nkl,wvt.Nkl,diagnostic)
                        testCase.verifyEqual(actual.spatialShape(:)',wvt.spatialMatrixSize,diagnostic)
                        testCase.verifyEqual(actual.spectralShape(:)',wvt.spectralMatrixSize,diagnostic)
                        testCase.verifyEqual(actual.j,wvt.j,diagnostic)
                        testCase.verifyEqual(actual.shouldAntialias,wvt.shouldAntialias,diagnostic)
                        testCase.verifyEqual(actual.kMode(:),wvt.kMode_wv(:),diagnostic)
                        testCase.verifyEqual(actual.lMode(:),wvt.lMode_wv(:),diagnostic)
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.k(:),wvt.k(:),1e-12,"k: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.l(:),wvt.l(:),1e-12,"l: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.Kh(:),wvt.Kh(:),1e-12,"Kh: " + diagnostic));
                        testCase.verifyEqual(uint64(actual.dftPrimaryIndices2D(:)),wvt.dftPrimaryIndices2D(:),diagnostic)
                        testCase.verifyEqual(uint64(actual.dftConjugateIndices2D(:)),wvt.dftConjugateIndices2D(:),diagnostic)
                        verifyHalfSpectrumMappings(testCase,actual,wvt,diagnostic)

                        deformationWavenumberSquared = double(j == 1)*wvt.f^2/(wvt.g*wvt.h);
                        testCase.verifyEqual(actual.coriolisFrequency,wvt.f,RelTol=1e-14)
                        testCase.verifyEqual(actual.deformationWavenumberSquared,deformationWavenumberSquared,RelTol=1e-14,AbsTol=eps)
                        factorPairs = { ...
                            "uFactor",wvt.UA0; ...
                            "vFactor",wvt.VA0};
                        for iFactor = 1:size(factorPairs,1)
                            name = factorPairs{iFactor,1};
                            value = actual.(name + "Real") + 1i*actual.(name + "Imag");
                            maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,value(:),factorPairs{iFactor,2}(:),1e-12,name + ": " + diagnostic));
                        end
                        realFactorPairs = { ...
                            "etaFactor",wvt.NA0; ...
                            "piFactor",wvt.PA0; ...
                            "psiFactor",wvt.A0_Psi_factor; ...
                            "qgpvFactor",wvt.A0_QGPV_factor; ...
                            "zetaZFactor",wvt.K2./(wvt.K2 + deformationWavenumberSquared); ...
                            "energyFactor",wvt.A0_TE_factor; ...
                            "enstrophyFactor",wvt.A0_TZ_factor};
                        for iFactor = 1:size(realFactorPairs,1)
                            name = realFactorPairs{iFactor,1};
                            expected = realFactorPairs{iFactor,2};
                            expected(wvt.Kh == 0) = 0;
                            maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.(name)(:),expected(:),1e-12,name + ": " + diagnostic));
                        end

                        index = (1:wvt.Nkl)';
                        A0 = 2e-5*sin(0.31*index) + 1i*1e-5*cos(0.17*(index+2));
                        A0 = reshape(A0,wvt.spectralMatrixSize);
                        A0(wvt.Kh == 0) = 0;
                        selfConjugate = wvt.dftPrimaryIndices2D == wvt.dftConjugateIndices2D;
                        A0(selfConjugate) = real(A0(selfConjugate));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.A0Real(:)+1i*actual.A0Imag(:),A0(:),1e-12,"A0 ordering: " + diagnostic));
                        wvt.A0 = A0;
                        [u,v,eta,pi,psi,qgpv,zetaZ,ssh] = wvt.variableWithName("u","v","eta","pi","psi","qgpv","zeta_z","ssh");
                        fieldPairs = { ...
                            "u",u; "v",v; "eta",eta; "pi",pi; ...
                            "psi",psi; "qgpv",qgpv; "zetaZ",zetaZ; "ssh",ssh};
                        for iField = 1:size(fieldPairs,1)
                            name = fieldPairs{iField,1};
                            maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.(name)(:),fieldPairs{iField,2}(:),1e-12,name + ": " + diagnostic));
                        end

                        expectedQGPVDerivatives = cat(3,qgpv,wvt.diffX(qgpv),wvt.diffY(qgpv));
                        expectedPsiDerivatives = cat(3,psi,wvt.diffX(psi),wvt.diffY(psi));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.qgpvWithDerivatives(:),expectedQGPVDerivatives(:),1e-12,"QGPV derivatives: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.psiWithDerivatives(:),expectedPsiDerivatives(:),1e-12,"psi derivatives: " + diagnostic));

                        projectedA0 = wvt.transformQGPVToWaveVortex(qgpv);
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.projectedA0Real(:)+1i*actual.projectedA0Imag(:),projectedA0(:),1e-12,"QGPV projection: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.evolvedA0Real(:)+1i*actual.evolvedA0Imag(:),A0(:),1e-12,"linear evolution: " + diagnostic));
                        nonlinearF0 = wvt.nonlinearFlux();
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.nonlinearF0Real(:)+1i*actual.nonlinearF0Imag(:),nonlinearF0(:),1e-12,"nonlinear PV tendency: " + diagnostic));

                        forcing = barotropicForcingDump(testCase.ForcingDump,definition);
                        testCase.verifyEqual(forcing.Nkl,wvt.Nkl,diagnostic)
                        forcingPairs = { ...
                            "nonlinear",forcingFlux(wvt,WVNonlinearAdvection(wvt)); ...
                            "damping",forcingFlux(wvt,WVAdaptiveDamping(wvt)); ...
                            "linear",forcingFlux(wvt,WVBottomFrictionLinear(wvt,r=2.5e-7)); ...
                            "quadratic",forcingFlux(wvt,WVBottomFrictionQuadratic(wvt,Cd=1.5e-3)); ...
                            "beta",forcingFlux(wvt,WVBetaPlanePVAdvection(wvt)); ...
                            "nonlinearDamping",forcingFlux(wvt,[WVNonlinearAdvection(wvt) WVAdaptiveDamping(wvt)]); ...
                            "nonlinearLinear",forcingFlux(wvt,[WVNonlinearAdvection(wvt) WVBottomFrictionLinear(wvt,r=2.5e-7)]); ...
                            "nonlinearQuadratic",forcingFlux(wvt,[WVNonlinearAdvection(wvt) WVBottomFrictionQuadratic(wvt,Cd=1.5e-3)]); ...
                            "betaDamping",forcingFlux(wvt,[WVBetaPlanePVAdvection(wvt) WVAdaptiveDamping(wvt)]); ...
                            "all",forcingFlux(wvt,[WVNonlinearAdvection(wvt) WVBottomFrictionQuadratic(wvt,Cd=1.5e-3) WVBottomFrictionLinear(wvt,r=2.5e-7) WVBetaPlanePVAdvection(wvt) WVAdaptiveDamping(wvt)])};
                        for iForcing = 1:size(forcingPairs,1)
                            name = forcingPairs{iForcing,1};
                            compiled = forcing.(name + "Real") + 1i*forcing.(name + "Imag");
                            maximumForcingRelativeError = max(maximumForcingRelativeError,verifyRelative(testCase,compiled(:),forcingPairs{iForcing,2}(:),1e-12,name + ": " + diagnostic));
                        end
                        fixedFlux = forcingFlux(wvt,[WVNonlinearAdvection(wvt) WVFixedAmplitudeForcing(wvt,name="fixed",A0_indices=uint64([2;4]),A0bar=[3e-6-2e-6i;-4e-6])]);
                        maximumForcingRelativeError = max(maximumForcingRelativeError,verifyRelative(testCase,forcing.fixedReal(:)+1i*forcing.fixedImag(:),fixedFlux(:),1e-12,"fixed amplitude: " + diagnostic));
                        testCase.verifyEqual(forcing.narrowReal,forcing.fixedReal,diagnostic)
                        testCase.verifyEqual(forcing.narrowImag,forcing.fixedImag,diagnostic)
                        testCase.verifyEqual(forcing.fieldReconstructionCount,1,diagnostic)
                        testCase.verifyEqual(forcing.fieldReuseCount,4,diagnostic)
                        testCase.verifyEqual(forcing.projectionCount,4,diagnostic)
                        testCase.verifyEqual(forcing.forcingCallCount,5,diagnostic)
                        testCase.verifyEqual(forcing.workspaceCapacityBytes,0,diagnostic)
                        if isequal(definition.Nxy,[9 6]) && j == 1 && shouldAntialias
                            endpointPairs = { ...
                                "rk4Endpoint","fixed"; ...
                                "rk23Endpoint","ode23"; ...
                                "rk45Endpoint","ode45"; ...
                                "rk78Endpoint","ode78"};
                            for iEndpoint = 1:size(endpointPairs,1)
                                name = endpointPairs{iEndpoint,1};
                                expectedEndpoint = forcingEndpoint(definition,A0,endpointPairs{iEndpoint,2});
                                compiledEndpoint = forcing.(name + "Real") + 1i*forcing.(name + "Imag");
                                maximumIntegrationRelativeError = max(maximumIntegrationRelativeError,verifyRelative(testCase,compiledEndpoint(:),expectedEndpoint(:),1e-12,name + ": " + diagnostic));
                            end
                        end

                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.totalEnergy,wvt.totalEnergy,1e-12,"spectral energy: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.totalEnergySpatiallyIntegrated,wvt.totalEnergySpatiallyIntegrated,1e-12,"spatial energy: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.totalEnstrophy,wvt.totalEnstrophy(),1e-12,"spectral enstrophy: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.totalEnstrophySpatiallyIntegrated,wvt.totalEnstrophySpatiallyIntegrated(),1e-12,"spatial enstrophy: " + diagnostic));
                        maximumRelativeError = max(maximumRelativeError,verifyRelative(testCase,actual.uvMax,wvt.uvMax,1e-12,"uvMax: " + diagnostic));

                        halfRows = (floor(wvt.Nx/2)+1)*wvt.Ny;
                        expectedScratchBytes = 4*halfRows*16 + 5*wvt.Nx*wvt.Ny*8;
                        testCase.verifyEqual(actual.scratchBytes,expectedScratchBytes,diagnostic)
                        testCase.verifyEqual(actual.persistentFullHermitianBytes,0,diagnostic)
                        testCase.verifyEqual(actual.planCount,3,diagnostic)
                        testCase.verifyEqual(string(actual.providerIdentifier),"reference-direct",diagnostic)
                        testCase.verifyEqual(string(actual.coefficientOrderingIdentifier),"matlab-kl-radial-k-l",diagnostic)
                        testCase.verifyEqual(string(actual.normalizationIdentifier),"A0-is-qgpv",diagnostic)
                        testCase.verifyEqual(string(actual.antialiasImplementationIdentifier),"transform-level-radial-two-thirds",diagnostic)
                    end
                end
            end
            fprintf("Barotropic QG MATLAB/C++ maximum relative error: %.3e\n",maximumRelativeError)
            fprintf("Barotropic QG MATLAB/C++ tolerance-vector maximum relative error: %.3e\n",maximumToleranceRelativeError)
            fprintf("Barotropic QG forcing MATLAB/C++ maximum relative error: %.3e\n",maximumForcingRelativeError)
            fprintf("Barotropic QG forcing integration MATLAB/C++ maximum relative error: %.3e\n",maximumIntegrationRelativeError)
        end
    end
end

function verifyHalfSpectrumMappings(testCase,actual,wvt,diagnostic)
layout = WVFourierStorageLayout(wvt,"hermitian-half",compressedDimension=1);
testCase.verifyEqual(uint64(actual.halfDirectRows(:))+1,layout.fourierRowsForDirectWVIndices(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfDirectWVIndices(:))+1,layout.directWVIndices(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfConjugatedRows(:))+1,layout.fourierRowsForConjugatedWVIndices(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfConjugatedWVIndices(:))+1,layout.conjugatedWVIndices(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfCompletionRows(:))+1,layout.hermitianCompletionRows(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfCompletionSourceRows(:))+1,layout.hermitianSourceRows(:),diagnostic)
testCase.verifyEqual(uint64(actual.halfSelfConjugateRows(:))+1,layout.selfConjugateFourierRows(:),diagnostic)
end

function relativeError = verifyElementwiseRelative(testCase,actual,expected,tolerance,diagnostic)
testCase.verifyEqual(size(actual),size(expected),diagnostic)
scale = max(abs(expected),realmin("double"));
relativeError = max(abs(actual-expected)./scale,[],"all");
testCase.verifyLessThanOrEqual(relativeError,tolerance,diagnostic)
end

function relativeError = verifyRelative(testCase,actual,expected,tolerance,diagnostic)
maximumExpected = max(abs(expected(:)),[],"all");
maximumError = max(abs(actual(:)-expected(:)),[],"all");
if maximumExpected == 0
    relativeError = 0;
    testCase.verifyLessThanOrEqual(maximumError,1e-14,diagnostic)
else
    relativeError = maximumError/maximumExpected;
    testCase.verifyLessThanOrEqual(relativeError,tolerance,diagnostic)
end
end

function actual = barotropicFixtureDump(definition)
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
executable = fullfile(repositoryRoot,"tools","compiled-kernel","build","WVBarotropicQGFixtureDump");
arguments = [definition.Nxy definition.Lxy definition.h definition.j definition.g definition.rotationRate definition.latitude definition.shouldAntialias definition.planetaryRadius];
command = sprintf('"%s" %s',executable,strjoin(compose("%.17g",arguments)," "));
[status,output] = systemWithoutMatlabRuntime(command);
if status ~= 0
    error("WaveVortexModel:BarotropicQGFixtureDumpFailed","%s",output);
end
actual = jsondecode(output);
end

function actual = barotropicToleranceDump(executable,definition,absoluteToleranceScale)
commandArguments = [definition.Nxy definition.Lxy definition.h definition.j definition.g definition.rotationRate definition.latitude definition.shouldAntialias definition.planetaryRadius absoluteToleranceScale];
command = shellQuote(executable) + " " + strjoin(compose("%.17g",commandArguments)," ");
[status,output] = systemWithoutMatlabRuntime(command);
if status ~= 0
    error("WaveVortexModel:BarotropicQGToleranceDumpFailed","%s",output);
end
actual = jsondecode(output);
end

function actual = barotropicForcingDump(executable,definition)
commandArguments = [definition.Nxy definition.Lxy definition.h definition.j definition.g definition.rotationRate definition.latitude definition.shouldAntialias definition.planetaryRadius];
command = shellQuote(executable) + " " + strjoin(compose("%.17g",commandArguments)," ");
[status,output] = systemWithoutMatlabRuntime(command);
if status ~= 0
    error("WaveVortexModel:BarotropicQGForcingDumpFailed","%s",output);
end
actual = jsondecode(output);
end

function F0 = forcingFlux(wvt,forcing)
wvt.removeAllForcing();
wvt.addForcing(forcing);
F0 = wvt.nonlinearFlux();
end

function A0 = forcingEndpoint(definition,initialA0,integrator)
wvt = WVTransformBarotropicQG( ...
    definition.Lxy,definition.Nxy, ...
    h=definition.h,j=definition.j,g=definition.g, ...
    rotationRate=definition.rotationRate, ...
    latitude=definition.latitude, ...
    shouldAntialias=definition.shouldAntialias, ...
    planetaryRadius=definition.planetaryRadius);
wvt.A0 = initialA0;
wvt.A0([2 4]) = [3e-6-2e-6i;-4e-6];
wvt.removeAllForcing();
wvt.addForcing([ ...
    WVNonlinearAdvection(wvt) ...
    WVBottomFrictionQuadratic(wvt,Cd=1.5e-3) ...
    WVBottomFrictionLinear(wvt,r=2.5e-7) ...
    WVBetaPlanePVAdvection(wvt) ...
    WVAdaptiveDamping(wvt) ...
    WVFixedAmplitudeForcing(wvt,name="fixed",A0_indices=uint64([2;4]),A0bar=[3e-6-2e-6i;-4e-6])]);
model = WVModel(wvt);
if integrator == "fixed"
    model.setupIntegrator(integratorType="fixed",deltaT=0.005);
else
    model.setupIntegrator(integratorType="adaptive",integrator=str2func(integrator),absTolerance=1e-10,relTolerance=1e-8);
end
model.integrateToTime(0.01,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
A0 = wvt.A0;
end

function value = shellQuote(value)
value = "'" + replace(string(value),"'","'""'""'") + "'";
end

function [status,output] = systemWithoutMatlabRuntime(command)
if isunix && ~ismac
    command = "env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH " + string(command);
end
[status,output] = system(command);
end
