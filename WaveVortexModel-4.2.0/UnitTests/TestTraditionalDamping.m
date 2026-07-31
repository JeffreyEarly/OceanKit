classdef TestTraditionalDamping < matlab.unittest.TestCase
    properties
        tempFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            testCase.tempFolder = string(tempname);
            mkdir(testCase.tempFolder);
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryFolder(testCase)
            if isfolder(testCase.tempFolder)
                rmdir(testCase.tempFolder,"s");
            end
        end
    end

    methods (Test)
        function verticalDampingSupportsConstantStratification(testCase)
            for isHydrostatic = [true false]
                wvt = TestTraditionalDamping.constantTransform(isHydrostatic=isHydrostatic,shouldAntialias=false);
                damping = WVVerticalDamping(wvt,nu=2e-4,kappa=3e-6);
                testCase.verifyEqual(damping.dLnN2,0)

                wvt.initWithRandomFlow(uvMax=0.01);
                if isHydrostatic
                    [Fu,Fv,Feta] = damping.addHydrostaticSpatialForcing(wvt,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize));
                    TestTraditionalDamping.verifyFinite(testCase,Fu,Fv,Feta);
                else
                    [Fu,Fv,Fw,Feta] = damping.addNonhydrostaticSpatialForcing(wvt,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize));
                    TestTraditionalDamping.verifyFinite(testCase,Fu,Fv,Fw,Feta);
                end
            end
        end

        function variableStratificationCorrectionIsPreserved(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*exp(z/1000);
            wvt = WVTransformHydrostatic([1e4,1e4,1000],[8,8,7],N2=N2,latitude=30,shouldAntialias=false);
            damping = WVVerticalDamping(wvt,nu=2e-4,kappa=3e-6);

            testCase.verifyEqual(damping.dLnN2,shiftdim(wvt.dLnN2,-2))
        end

        function directResolutionConversionPreservesCoefficients(testCase)
            source = TestTraditionalDamping.constantTransform(shouldAntialias=false);
            target = source.waveVortexTransformWithResolution([12,12,7]);
            horizontal = WVHorizontalDamping(source,nu=0.125,kappa=2.5e-5);
            vertical = WVVerticalDamping(source,nu=0.375,kappa=7.5e-5);

            horizontalConverted = horizontal.forcingWithResolutionOfTransform(target);
            verticalConverted = vertical.forcingWithResolutionOfTransform(target);

            TestTraditionalDamping.verifyCoefficients(testCase,horizontalConverted,horizontal.nu,horizontal.kappa)
            TestTraditionalDamping.verifyCoefficients(testCase,verticalConverted,vertical.nu,vertical.kappa)
            testCase.verifyTrue(horizontalConverted.isClosure)
            testCase.verifyTrue(verticalConverted.isClosure)
        end

        function transformResolutionConversionPreservesCoefficients(testCase)
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=false);
            TestTraditionalDamping.addClosures(wvt);

            converted = wvt.waveVortexTransformWithResolution([12,12,7]);

            TestTraditionalDamping.verifyTransformClosures(testCase,converted)
        end

        function explicitAntialiasingPreservesCoefficients(testCase)
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=true);
            TestTraditionalDamping.addClosures(wvt);

            converted = wvt.waveVortexTransformWithExplicitAntialiasing();

            testCase.verifyFalse(converted.shouldAntialias)
            testCase.verifyTrue(converted.hasForcingWithName("antialias filter"))
            TestTraditionalDamping.verifyTransformClosures(testCase,converted)
        end

        function horizontalTendencyMatchesResolvedLaplacian(testCase)
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=false);
            nu = 0.125;
            kappa = 2.5e-5;
            damping = WVHorizontalDamping(wvt,nu=nu,kappa=kappa);
            k = 2*pi/wvt.Lx;
            l = 2*pi/wvt.Ly;
            U = sin(k*wvt.X);
            V = cos(l*wvt.Y);
            Eta = sin(k*wvt.X+l*wvt.Y);
            wvt.initWithUVEta(U,V,Eta);

            [Fu,Fv,Fw,Feta] = damping.addNonhydrostaticSpatialForcing(wvt,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize));

            testCase.verifyEqual(Fu,-nu*k*k*wvt.u,AbsTol=1e-12)
            testCase.verifyEqual(Fv,-nu*l*l*wvt.v,AbsTol=1e-12)
            testCase.verifyEqual(Fw,nu*(wvt.diffX(wvt.w,n=2)+wvt.diffY(wvt.w,n=2)),AbsTol=1e-12)
            testCase.verifyEqual(Feta,-kappa*(k*k+l*l)*wvt.eta,AbsTol=1e-12)
        end

        function verticalTendencyMatchesResolvedLaplacian(testCase)
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=false);
            nu = 0.375;
            kappa = 7.5e-5;
            damping = WVVerticalDamping(wvt,nu=nu,kappa=kappa);
            j = 1;
            m = j*pi/wvt.Lz;
            wvt.initWithWaveModes(kMode=1,lMode=0,j=j,phi=0,u=0.01,sign=1);

            [Fu,Fv,Fw,Feta] = damping.addNonhydrostaticSpatialForcing(wvt,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize));

            testCase.verifyEqual(Fu,-nu*m*m*wvt.u,AbsTol=1e-12)
            testCase.verifyEqual(Fv,-nu*m*m*wvt.v,AbsTol=1e-12)
            testCase.verifyEqual(Fw,-nu*m*m*wvt.w,AbsTol=1e-12)
            testCase.verifyEqual(Feta,-kappa*m*m*wvt.eta,AbsTol=1e-12)
        end

        function netCDFRoundTripPreservesClosures(testCase)
            path = fullfile(testCase.tempFolder,"traditional-damping.nc");
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=true);
            TestTraditionalDamping.addClosures(wvt);
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            ncfile.close();

            restored = WVTransform.waveVortexTransformFromFile(char(path));

            TestTraditionalDamping.verifyTransformClosures(testCase,restored)
            testCase.verifyEqual(restored.forcingWithName("vertical scalar diffusivity").dLnN2,0)
        end

        function modelRestartPreservesClosures(testCase)
            path = fullfile(testCase.tempFolder,"traditional-damping-model.nc");
            wvt = TestTraditionalDamping.constantTransform(shouldAntialias=true);
            TestTraditionalDamping.addClosures(wvt);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            restoredModel = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()restoredModel.closeNetCDFFile());
            TestTraditionalDamping.verifyTransformClosures(testCase,restoredModel.wvt)
            vertical = restoredModel.wvt.forcingWithName("vertical scalar diffusivity");
            [Fu,Fv,Fw,Feta] = vertical.addNonhydrostaticSpatialForcing(restoredModel.wvt,zeros(restoredModel.wvt.spatialMatrixSize),zeros(restoredModel.wvt.spatialMatrixSize),zeros(restoredModel.wvt.spatialMatrixSize),zeros(restoredModel.wvt.spatialMatrixSize));
            TestTraditionalDamping.verifyFinite(testCase,Fu,Fv,Fw,Feta);
            restoredModel.closeNetCDFFile();
            clear cleanup
        end
    end

    methods (Static, Access=private)
        function wvt = constantTransform(options)
            arguments
                options.isHydrostatic (1,1) logical = false
                options.shouldAntialias (1,1) logical = false
            end
            wvt = WVTransformConstantStratification([1e4,1e4,1000],[8,8,5],N0=5.2e-3,latitude=30,isHydrostatic=options.isHydrostatic,shouldAntialias=options.shouldAntialias);
        end

        function addClosures(wvt)
            wvt.addForcing(WVHorizontalDamping(wvt,nu=0.125,kappa=2.5e-5));
            wvt.addForcing(WVVerticalDamping(wvt,nu=0.375,kappa=7.5e-5));
        end

        function verifyTransformClosures(testCase,wvt)
            horizontal = wvt.forcingWithName("horizontal scalar diffusivity");
            vertical = wvt.forcingWithName("vertical scalar diffusivity");
            TestTraditionalDamping.verifyCoefficients(testCase,horizontal,0.125,2.5e-5)
            TestTraditionalDamping.verifyCoefficients(testCase,vertical,0.375,7.5e-5)
            testCase.verifyTrue(horizontal.isClosure)
            testCase.verifyTrue(vertical.isClosure)
        end

        function verifyCoefficients(testCase,forcing,nu,kappa)
            testCase.verifyEqual(forcing.nu,nu)
            testCase.verifyEqual(forcing.kappa,kappa)
        end

        function verifyFinite(testCase,varargin)
            for iVariable = 1:length(varargin)
                testCase.verifyTrue(all(isfinite(varargin{iVariable}),"all"))
            end
        end
    end
end
