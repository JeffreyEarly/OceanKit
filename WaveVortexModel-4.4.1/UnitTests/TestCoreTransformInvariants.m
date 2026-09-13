classdef TestCoreTransformInvariants < matlab.unittest.TestCase
    properties
        wvt
        transformType
        gridType
    end

    properties (ClassSetupParameter)
        transform = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG", ...
            barotropicQG="barotropicQG")
        grid = struct(even=[8 6 9],odd=[9 7 8])
    end

    methods (TestClassSetup)
        function createTransform(testCase,transform,grid)
            testCase.transformType = string(transform);
            testCase.gridType = string(matlab.lang.makeValidName(sprintf("%dx%dx%d",grid)));
            testCase.wvt = testCase.transformForConfiguration(transform,grid,shouldAntialias=false);
        end
    end

    methods (Test, TestTags = "full")
        function geometryConjugacyAndIndexBijections(testCase)
            wvt = testCase.wvt;
            diagnostic = testCase.diagnostic("geometry, conjugacy, and indexing");

            degreesOfFreedom = WVGeometryDoublyPeriodic.degreesOfFreedomForRealMatrix(wvt.Nx,wvt.Ny,conjugateDimension=wvt.conjugateDimension);
            testCase.verifyEqual(sum(degreesOfFreedom,"all"),wvt.Nx*wvt.Ny,diagnostic)

            nyquistMask = WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny);
            if mod(wvt.Nx,2) == 0
                testCase.verifyTrue(all(nyquistMask(wvt.Nx/2+1,:),"all"),diagnostic)
            else
                testCase.verifyFalse(any(nyquistMask,"all") && mod(wvt.Ny,2) == 1,diagnostic)
                testCase.verifyTrue(any(wvt.kMode_wv == floor(wvt.Nx/2)),diagnostic)
            end
            if mod(wvt.Ny,2) == 0
                testCase.verifyTrue(all(nyquistMask(:,wvt.Ny/2+1),"all"),diagnostic)
            else
                testCase.verifyTrue(any(wvt.lMode_wv == floor(wvt.Ny/2)),diagnostic)
            end

            spectralIndices = (1:prod(wvt.spectralMatrixSize))';
            [kMode,lMode,jMode] = wvt.modeNumberFromIndex(spectralIndices);
            testCase.verifyEqual(wvt.indexFromModeNumber(kMode,lMode,jMode),spectralIndices,diagnostic)
            conjugateCandidate = find((kMode ~= 0 | lMode ~= 0) & wvt.isValidConjugateModeNumber(-kMode,-lMode,jMode),1);
            testCase.assertNotEmpty(conjugateCandidate,diagnostic)
            testCase.verifyEqual(wvt.indexFromModeNumber(-kMode(conjugateCandidate),-lMode(conjugateCandidate),jMode(conjugateCandidate)),spectralIndices(conjugateCandidate),diagnostic)
            reordered = spectralIndices([numel(spectralIndices); 1; min(3,numel(spectralIndices))]);
            [kReordered,lReordered,jReordered] = wvt.modeNumberFromIndex(reordered);
            testCase.verifyEqual(wvt.indexFromModeNumber(kReordered,lReordered,jReordered),reordered,diagnostic)
            testCase.verifyFalse(wvt.isValidPrimaryModeNumber(kMode(1),lMode(1),max(wvt.j)+1),diagnostic)
            testCase.verifyError(@()wvt.modeNumberFromIndex(prod(wvt.spectralMatrixSize)+1),"MATLAB:validators:mustBeLessThanOrEqual")
            testCase.verifyError(@()wvt.klModeNumberFromIndex(wvt.Nkl+1),"MATLAB:validators:mustBeLessThanOrEqual")

            seedRandomNumberGenerator(testCase,14002);
            spatialField = randn(wvt.spatialMatrixSize);
            dftField = wvt.transformFromSpatialDomainToDFTGrid(spatialField);
            testCase.verifyTrue(WVGeometryDoublyPeriodic.isHermitian(dftField),diagnostic)
            testCase.verifyLessThanOrEqual(abs(imag(dftField(1,1))),10*eps,diagnostic)
            dftField(1,1) = dftField(1,1)+1i;
            testCase.verifyFalse(WVGeometryDoublyPeriodic.isHermitian(dftField),diagnostic)
        end

        function physicalSpectralRoundTripsAndQuadraticInvariants(testCase)
            wvt = testCase.wvt;
            diagnostic = testCase.diagnostic("round trip, energy, and enstrophy");
            seedRandomNumberGenerator(testCase,14003);
            wvt.initWithRandomFlow(uvMax=0.01);

            switch testCase.transformType
                case {"constantHydrostatic","hydrostatic"}
                    [u,v,eta] = wvt.variableWithName("u","v","eta");
                    [Ap,Am,A0] = wvt.transformUVEtaToWaveVortex(u,v,eta);
                    testCase.verifyCoefficients(Ap,Am,A0,5e-10,diagnostic)
                case {"constantNonhydrostatic","boussinesq"}
                    [u,v,w,eta] = wvt.variableWithName("u","v","w","eta");
                    [Ap,Am,A0] = wvt.transformUVWEtaToWaveVortex(u,v,w,eta);
                    coefficientTolerance = 5e-10;
                    if testCase.transformType == "boussinesq"
                        coefficientTolerance = 1e-3;
                    end
                    testCase.verifyCoefficients(Ap,Am,A0,coefficientTolerance,diagnostic)
                otherwise
                    qgpv = wvt.variableWithName("qgpv");
                    A0 = wvt.transformQGPVToWaveVortex(qgpv);
                    testCase.verifyRelative(A0,wvt.A0,5e-10,diagnostic)
            end

            testCase.verifyRelative(wvt.totalEnergySpatiallyIntegrated,wvt.totalEnergy,1e-3,diagnostic)
            testCase.verifyRelative(wvt.totalEnstrophySpatiallyIntegrated,wvt.totalEnstrophy,1e-3,diagnostic)
        end

        function analyticalHorizontalDifferentiation(testCase)
            wvt = testCase.wvt;
            diagnostic = testCase.diagnostic("horizontal differentiation");
            X = wvt.X;
            Y = wvt.Y;
            k = 2*pi*min(2,floor((wvt.Nx-1)/2))/wvt.Lx;
            l = 2*pi*min(2,floor((wvt.Ny-1)/2))/wvt.Ly;
            field = cos(k*X+l*Y);
            testCase.verifyRelative(wvt.diffX(field),-k*sin(k*X+l*Y).*ones(size(field)),5e-12,diagnostic)
            testCase.verifyRelative(wvt.diffY(field),-l*sin(k*X+l*Y).*ones(size(field)),5e-12,diagnostic)
            testCase.verifyRelative(wvt.diffX(field,n=4),k^4*field,5e-12,diagnostic)
            testCase.verifyRelative(wvt.diffY(field,n=4),l^4*field,5e-12,diagnostic)
        end

        function resolutionAndAntialiasConversionsPreserveCommonModes(testCase)
            wvt = testCase.wvt;
            diagnostic = testCase.diagnostic("resolution conversion");
            seedRandomNumberGenerator(testCase,14004);
            wvt.initWithRandomFlow(uvMax=0.01);
            wvt.t0 = 11;
            wvt.t = 13;

            if numel(wvt.spatialMatrixSize) == 3
                largerSize = [wvt.Nx+2 wvt.Ny+2 wvt.Nz+1];
                smallerSize = [max(5,wvt.Nx-2) max(5,wvt.Ny-2) max(5,wvt.Nz-1)];
            else
                largerSize = [wvt.Nx+2 wvt.Ny+2];
                smallerSize = [max(5,wvt.Nx-2) max(5,wvt.Ny-2)];
            end
            larger = wvt.waveVortexTransformWithResolution(largerSize);
            smaller = wvt.waveVortexTransformWithResolution(smallerSize);
            testCase.verifyConversion(wvt,larger,diagnostic)
            testCase.verifyConversion(wvt,smaller,diagnostic)

            if testCase.gridType == "x8x6x9" && any(testCase.transformType == ["constantHydrostatic" "constantNonhydrostatic" "hydrostatic" "boussinesq"])
                antialiased = testCase.transformForConfiguration(testCase.transformType,[8 6 9],shouldAntialias=true);
                seedRandomNumberGenerator(testCase,14005);
                antialiased.initWithRandomFlow(uvMax=0.01);
                explicit = antialiased.waveVortexTransformWithExplicitAntialiasing();
                testCase.verifyFalse(explicit.shouldAntialias,diagnostic)
                testCase.verifyTrue(explicit.hasForcingWithName("antialias filter"),diagnostic)
                testCase.verifyConversion(antialiased,explicit,diagnostic)
            end
        end
    end

    methods (Access = private)
        function verifyCoefficients(testCase,Ap,Am,A0,tolerance,diagnostic)
            testCase.verifyRelative(Ap,testCase.wvt.Ap,tolerance,diagnostic)
            testCase.verifyRelative(Am,testCase.wvt.Am,tolerance,diagnostic)
            testCase.verifyRelative(A0,testCase.wvt.A0,tolerance,diagnostic)
        end

        function verifyConversion(testCase,source,target,diagnostic)
            testCase.verifyEqual(class(target),class(source),diagnostic)
            testCase.verifyEqual([target.Lx target.Ly target.Lz],[source.Lx source.Ly source.Lz],diagnostic)
            testCase.verifyEqual(target.t0,source.t0,diagnostic)
            testCase.verifyEqual(target.t,source.t,diagnostic)
            if isa(source,"WVGeometryDoublyPeriodicBarotropic")
                testCase.verifyEqual(target.j,source.j,diagnostic)
            end

            sourceK = repmat(source.kMode_wv.',source.Nj,1);
            sourceL = repmat(source.lMode_wv.',source.Nj,1);
            sourceJ = repmat(source.j,1,source.Nkl);
            targetK = repmat(target.kMode_wv.',target.Nj,1);
            targetL = repmat(target.lMode_wv.',target.Nj,1);
            targetJ = repmat(target.j,1,target.Nkl);
            [isCommon,sourceIndex] = ismember([targetK(:),targetL(:),targetJ(:)],[sourceK(:),sourceL(:),sourceJ(:)],"rows");
            for coefficientName = ["A0" "Ap" "Am"]
                if isempty(source.(coefficientName))
                    testCase.verifyEmpty(target.(coefficientName),diagnostic)
                else
                    targetCoefficient = target.(coefficientName);
                    sourceCoefficient = source.(coefficientName);
                    if numel(sourceCoefficient) ~= prod(source.spectralMatrixSize)
                        testCase.verifyEqual(targetCoefficient,sourceCoefficient,diagnostic)
                        continue
                    end
                    targetCoefficient = targetCoefficient(:);
                    sourceCoefficient = sourceCoefficient(:);
                    testCase.verifyEqual(targetCoefficient(isCommon),sourceCoefficient(sourceIndex(isCommon)),diagnostic)
                    testCase.verifyEqual(targetCoefficient(~isCommon),zeros(size(targetCoefficient(~isCommon))),diagnostic)
                end
            end
        end

        function text = diagnostic(testCase,invariant)
            text = sprintf("%s on %s grid: %s",testCase.transformType,testCase.gridType,invariant);
        end

        function verifyRelative(testCase,actual,expected,tolerance,diagnostic)
            relativeError = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
            testCase.verifyLessThanOrEqual(relativeError,tolerance,diagnostic)
        end
    end

    methods (Static, Access = private)
        function wvt = transformForConfiguration(transform,Nxyz,options)
            arguments
                transform string
                Nxyz (1,3) double
                options.shouldAntialias (1,1) logical
            end
            Lxyz = [6e3 4e3 1e3];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transform
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=options.shouldAntialias);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=options.shouldAntialias);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=options.shouldAntialias);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=options.shouldAntialias);
                case "stratifiedQG"
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=options.shouldAntialias);
                case "barotropicQG"
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=options.shouldAntialias,j=0);
            end
        end
    end
end
