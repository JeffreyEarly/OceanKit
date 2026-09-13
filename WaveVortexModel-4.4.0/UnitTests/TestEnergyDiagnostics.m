classdef TestEnergyDiagnostics < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function qgEnergyUsesOnlyPresentCoefficientFamilies(testCase)
            for family = ["barotropic","stratified-qg"]
                wvt = TestEnergyDiagnostics.transform(family);
                wvt.addOperation(wvt.operationForKnownVariable('energy'));
                testCase.verifyFalse(wvt.totalFlowComponent.hasWaveComponent);
                beforeWaves = {wvt.Ap,wvt.Am};
                testCase.verifyEqual(wvt.energy,0);
                TestEnergyDiagnostics.populate(wvt);
                testCase.verifyEqual({wvt.Ap,wvt.Am},beforeWaves);
                testCase.verifySize(wvt.energy,[1 1]);
                testCase.verifyEqual(wvt.energy,wvt.totalEnergy);
                testCase.verifyGreaterThan(wvt.energy,0);
                for componentName = wvt.flowComponentNames
                    component = wvt.flowComponentWithName(componentName);
                    operation = wvt.operationForKnownVariable('energy',flowComponent=component);
                    testCase.verifyEqual(operation.compute(wvt),wvt.totalEnergyOfFlowComponent(component),RelTol=1e-14);
                end
                mask = wvt.geostrophicComponent.maskA0 & wvt.Kh>wvt.dk;
                component = WVFlowComponent(wvt,maskA0=mask);
                component.name = 'selected geostrophic modes';
                component.abbreviatedName = 'selected';
                operation = wvt.operationForKnownVariable('energy',flowComponent=component);
                expected = sum(wvt.A0_TE_factor(:).*mask(:).*abs(wvt.A0(:)).^2);
                testCase.verifyEqual(operation.compute(wvt),expected,RelTol=1e-14);
                % A requested family absent from the transform contributes
                % zero, even when a custom component selects that family.
                component.maskAp = 1; component.maskAm = 1;
                operation = wvt.operationForKnownVariable('energy',flowComponent=component);
                testCase.verifyEqual(operation.compute(wvt),expected,RelTol=1e-14);
                empty = WVFlowComponent(wvt);
                empty.name = 'empty'; empty.abbreviatedName = 'empty';
                operation = wvt.operationForKnownVariable('energy',flowComponent=empty);
                testCase.verifyEqual(operation.compute(wvt),0);
            end
        end

        function observerPreservesRegisteredEnergyOperation(testCase)
            wvt = TestEnergyDiagnostics.transform("barotropic");
            annotation = WVVariableAnnotation('energy',{},'m3 s-2','user energy');
            operation = WVOperation('energy',annotation,@(~)42);
            wvt.addOperation(operation);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            testCase.verifyWarningFree(@()WVEulerianFields(model,fieldNames="energy"));
            testCase.verifyEqual(wvt.energy,42);
            testCase.verifyTrue(wvt.propertyAnnotationWithName('energy').modelOp==operation);
        end

        function waveEnergyRetainsExistingSummationAndMasks(testCase)
            for family = ["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                wvt = TestEnergyDiagnostics.transform(family);
                TestEnergyDiagnostics.populate(wvt);
                components = [wvt.totalFlowComponent,wvt.geostrophicComponent,wvt.waveComponent,wvt.inertialComponent,wvt.mdaComponent];
                for component = components
                    [~,~,mask] = WVTransform.optimizedTransformsForFlowComponent(wvt.totalFlowComponent,component);
                    expected = sum(wvt.Apm_TE_factor(:).*(mask.Ap(:).*abs(wvt.Ap(:)).^2+mask.Am(:).*abs(wvt.Am(:)).^2) ...
                        + wvt.A0_TE_factor(:).*(mask.A0(:).*abs(wvt.A0(:)).^2));
                    operation = wvt.operationForKnownVariable('energy',flowComponent=component);
                    testCase.verifyEqual(operation.compute(wvt),expected);
                end
            end
        end

        function energyOutputSurvivesMatlabReload(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for isLinear = [true false]
                for family = ["barotropic","stratified-qg","constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                    wvt = TestEnergyDiagnostics.transform(family);
                    TestEnergyDiagnostics.populate(wvt);
                    wvt.addOperation(wvt.operationForKnownVariable('energy'));
                    expected = wvt.energy;
                    wvt.setForcing(WVAntialiasing(wvt,Nj=max(wvt.j)+1));
                    model = WVModel(wvt,shouldUseLinearDynamics=isLinear);
                    model.eulerianObservingSystem.addNetCDFOutputVariables('energy');
                    path = fullfile(fixture.Folder,family+"-"+isLinear+".nc");
                    output = model.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true);
                    output.outputTimesForIntegrationPeriod(0,1);
                    output.writeTimeStepToOutputFile(0);
                    model.closeNetCDFFile();
                    testCase.verifyEqual(ncread(path,'/wave-vortex/energy'),expected);
                    testCase.verifyEqual(ncreadatt(path,'/wave-vortex/energy','units'),'m3 s-2');
                    restored = testCase.verifyWarningFree(@()WVModel.modelFromFile(path));
                    cleanup = onCleanup(@()restored.closeNetCDFFile());
                    testCase.verifyEqual(restored.wvt.energy,expected,RelTol=1e-14);
                    restored.setupIntegrator(integratorType="fixed",deltaT=.25);
                    restored.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                    restored.closeNetCDFFile();
                    clear cleanup
                    testCase.verifyEqual(ncread(path,'/wave-vortex/t'),[0;.5;1]);
                    % Preserve initial-only linear output and a time series in
                    % nonlinear mode; zero forcing holds coefficients constant.
                    if isLinear
                        testCase.verifyEqual(ncread(path,'/wave-vortex/energy'),expected,RelTol=1e-14);
                    else
                        testCase.verifyEqual(ncread(path,'/wave-vortex/energy'),repmat(expected,3,1),RelTol=1e-14);
                    end
                end
            end
        end
    end
    methods (Static,Access=private)
        function wvt = transform(family)
            switch family
                case "barotropic"
                    wvt = WVTransformBarotropicQG([17000 11000],[8 6],j=1,shouldAntialias=false);
                case "stratified-qg"
                    wvt = WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
                case {"constant-hydrostatic","constant-nonhydrostatic"}
                    wvt = WVTransformConstantStratification([17000 11000 1000],[8 6 9],N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=false);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic([17000 11000 1000],[8 6 9],Nj=4,N2=@(z)1e-4*exp(z/700),shouldAntialias=false);
                case "boussinesq"
                    wvt = WVTransformBoussinesq([17000 11000 1000],[8 6 9],Nj=4,N2=@(z)1e-4*exp(z/700),shouldAntialias=false);
            end
        end
        function populate(wvt)
            n = reshape(1:numel(wvt.A0),size(wvt.A0));
            wvt.A0 = 1e-6*complex(sin(.17*n),cos(.23*n)).*wvt.totalFlowComponent.maskA0;
            if wvt.totalFlowComponent.hasWaveComponent
                wvt.Ap = 1e-4*complex(sin(.13*n),cos(.19*n)).*wvt.totalFlowComponent.maskAp;
                wvt.Am = 1e-4*complex(sin(.29*n),cos(.31*n)).*wvt.totalFlowComponent.maskAm;
            end
        end
    end
end
