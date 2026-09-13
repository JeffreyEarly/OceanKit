classdef TestSQGTracerConvenience < matlab.unittest.TestCase
    properties
        tempFolder string
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.tempFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="smoke")
        function convenienceUsesTransformVelocityCapability(testCase)
            transforms = {
                WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false)
                WVTransformBarotropicQG([17000 11000],[8 6],shouldAntialias=false)
                WVTransformHydrostatic([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false)
                WVTransformBoussinesq([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false)
                };
            expectedXYOnly = [true true false false];
            for index = 1:numel(transforms)
                wvt = transforms{index};
                model = WVModel(wvt,shouldUseLinearDynamics=true);
                phi = reshape(1:prod(wvt.spatialMatrixSize),wvt.spatialMatrixSize);
                model.addTracer(phi,"dye");
                tracer = model.fluxedObservingSystemWithName("dye");
                testCase.verifyEqual(tracer.isXYOnly,expectedXYOnly(index))
                testCase.verifyEqual(tracer.isXYOnly,~wvt.hasVariableWithName('w'))
                testCase.verifyEqual(size(tracer.phi),size(phi))
                testCase.verifyEqual(tracer.phi,phi)
            end
        end

        function sqgConvenienceMatchesExplicitHorizontalRHSAndIntegration(testCase)
            convenience = TestSQGTracerConvenience.sqgModel("convenience");
            explicit = TestSQGTracerConvenience.sqgModel("explicit-horizontal");
            convenienceTracer = convenience.fluxedObservingSystemWithName("dye");
            explicitTracer = explicit.fluxedObservingSystemWithName("dye");

            testCase.verifyTrue(convenienceTracer.isXYOnly)
            testCase.verifyTrue(explicitTracer.isXYOnly)
            convenienceRHS = convenienceTracer.fluxAtTime(0,convenienceTracer.initialConditions());
            explicitRHS = explicitTracer.fluxAtTime(0,explicitTracer.initialConditions());
            testCase.verifyGreaterThan(max(abs(convenienceRHS{1}),[],"all"),0)
            testCase.verifyEqual(convenienceRHS,explicitRHS,AbsTol=1e-14)

            convenience.setupIntegrator(integratorType="fixed",deltaT=.25);
            explicit.setupIntegrator(integratorType="fixed",deltaT=.25);
            convenience.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            explicit.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.verifyEqual(convenience.tracer("dye"),explicit.tracer("dye"),AbsTol=1e-13)
            testCase.verifyEqual(convenience.wvt.A0,explicit.wvt.A0,AbsTol=1e-13)
        end

        function sqgConvenienceTracerRestartsWithXYZLayout(testCase)
            uninterrupted = TestSQGTracerConvenience.sqgModel("explicit-horizontal");
            checkpoint = TestSQGTracerConvenience.sqgModel("convenience");
            path = fullfile(testCase.tempFolder,"sqg-tracer-restart.nc");
            checkpoint.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true);
            uninterrupted.setupIntegrator(integratorType="fixed",deltaT=.25);
            checkpoint.setupIntegrator(integratorType="fixed",deltaT=.25);

            uninterrupted.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            checkpoint.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            checkpoint.closeNetCDFFile();

            info = ncinfo(path,"/wave-vortex/dye");
            testCase.verifyEqual(string({info.Dimensions.Name}),["/x" "/y" "/z" "t"])
            testCase.verifyEqual([info.Dimensions(1:3).Length],[8 6 9])

            resumed = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()resumed.closeNetCDFFile());
            restoredTracer = resumed.fluxedObservingSystemWithName("dye");
            testCase.verifyTrue(restoredTracer.isXYOnly)
            testCase.verifyEqual(size(restoredTracer.phi),[8 6 9])
            resumed.setupIntegrator(integratorType="fixed",deltaT=.25);
            resumed.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.verifyEqual(resumed.tracer("dye"),uninterrupted.tracer("dye"),AbsTol=1e-13)
            testCase.verifyEqual(resumed.wvt.A0,uninterrupted.wvt.A0,AbsTol=1e-13)
            resumed.closeNetCDFFile();
            clear cleanup
        end
    end

    methods (Static, Access=private)
        function model = sqgModel(tracerConstruction)
            wvt = WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
            index = reshape(1:numel(wvt.A0),size(wvt.A0));
            wvt.A0 = 1e-6*complex(sin(.17*index),cos(.23*index)).*(wvt.Kh>0);
            wvt.addForcing(WVAdaptiveDamping(wvt));
            model = WVModel(wvt);
            phi = .3+.2*sin(2*pi*wvt.X/wvt.Lx).*cos(2*pi*wvt.Y/wvt.Ly).*exp(wvt.Z/1000);
            if tracerConstruction == "convenience"
                model.addTracer(phi,"dye");
            else
                model.addFluxedObservingSystem(WVTracer(model,name="dye",phi=phi,isXYOnly=true));
            end
        end
    end
end
