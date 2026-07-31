classdef TestWVPseudoTopographicWaveGenerationProduction < matlab.unittest.TestCase
    % Verify variable stratification, resolution rebuilding, and restart behavior.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
        end
    end

    methods (Test)
        function variableStratificationMatchesOracleAndIdentities(testCase)
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(7142,"twister")
            for shouldAntialias = [false true]
                wvt = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],shouldAntialias);
                terrain = TestWVPseudoTopographicWaveGenerationProduction.topography(wvt);
                forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02+0.015i],frequency=1.405e-4);
                for t = [0 813]
                    wvt.t = t;
                    velocity = forcing.barotropicVelocityAtTime(t);
                    [directFp,directFm] = referenceBoundaryProjection(wvt,terrain,velocity);
                    [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                    testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(Fp,directFp),1e-10)
                    testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(Fm,directFm),1e-10)
                    testCase.verifyEqual(F0,zeros(size(wvt.A0)))
                end

                wvt.Ap = (randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
                wvt.Am = (randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
                wvt.A0(:) = 0;
                diagnostics = TestWVPseudoTopographicWaveGenerationProduction.sourceDiagnostics(wvt,forcing);
                testCase.verifyEqual(diagnostics.F0,zeros(size(wvt.A0)))
                testCase.verifyEqual(diagnostics.modalQGPVSource,zeros(size(wvt.A0)))
                testCase.verifyLessThanOrEqual(diagnostics.qgpvNorm,1e-10)
                testCase.verifyLessThanOrEqual(diagnostics.powerError,1e-10)
            end
            clear randomStateCleanup
        end

        function resolutionConversionRebuildsEquivalentForcing(testCase)
            sourceTransform = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],false);
            terrain = TestWVPseudoTopographicWaveGenerationProduction.topography(sourceTransform);
            forcing = WVPseudoTopographicWaveGeneration(sourceTransform,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02],frequency=1.31e-4,rampDuration=200,startTime=50,shouldAvoidAdaptiveDamping=false,maximumForcedHorizontalWavenumber=2.5e-4,maximumForcedVerticalMode=3,name="converted forcing");
            targetTransform = TestWVPseudoTopographicWaveGenerationProduction.createTransform([12 10 7],false);
            converted = forcing.forcingWithResolutionOfTransform(targetTransform);
            expectedTerrain = TestWVPseudoTopographicWaveGenerationProduction.topography(targetTransform);
            expected = WVPseudoTopographicWaveGeneration(targetTransform,topographicHeight=expectedTerrain,barotropicVelocityAmplitude=forcing.barotropicVelocityAmplitude,frequency=forcing.frequency,rampDuration=forcing.rampDuration,startTime=forcing.startTime,shouldAvoidAdaptiveDamping=forcing.shouldAvoidAdaptiveDamping,maximumForcedHorizontalWavenumber=forcing.maximumForcedHorizontalWavenumber,maximumForcedVerticalMode=forcing.maximumForcedVerticalMode,name=string(forcing.name));

            testCase.verifyEqual(converted.topographicHeight,expectedTerrain,AbsTol=1e-12)
            testCase.verifyEqual(converted.barotropicVelocityAmplitude,forcing.barotropicVelocityAmplitude)
            testCase.verifyEqual(converted.frequency,forcing.frequency)
            testCase.verifyEqual(converted.darwinSymbol,forcing.darwinSymbol)
            testCase.verifyEqual(converted.rampDuration,forcing.rampDuration)
            testCase.verifyEqual(converted.startTime,forcing.startTime)
            testCase.verifyEqual(converted.shouldAvoidAdaptiveDamping,forcing.shouldAvoidAdaptiveDamping)
            testCase.verifyEqual(converted.maximumForcedHorizontalWavenumber,forcing.maximumForcedHorizontalWavenumber)
            testCase.verifyEqual(converted.maximumForcedVerticalMode,forcing.maximumForcedVerticalMode)
            testCase.verifyEqual(string(converted.name),string(forcing.name))
            targetTransform.t = 467;
            [convertedFp,convertedFm,convertedF0] = converted.addSpectralForcing(targetTransform,zeros(size(targetTransform.Ap)),zeros(size(targetTransform.Am)),zeros(size(targetTransform.A0)));
            [expectedFp,expectedFm,expectedF0] = expected.addSpectralForcing(targetTransform,zeros(size(targetTransform.Ap)),zeros(size(targetTransform.Am)),zeros(size(targetTransform.A0)));
            testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(convertedFp,expectedFp),1e-12)
            testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(convertedFm,expectedFm),1e-12)
            testCase.verifyEqual(convertedF0,expectedF0)

            [xHigh,~] = ndgrid(targetTransform.x,targetTransform.y);
            highTerrain = expectedTerrain+7*cos(2*pi*5*xHigh/targetTransform.Lx);
            highForcing = WVPseudoTopographicWaveGeneration(targetTransform,topographicHeight=highTerrain,barotropicVelocityAmplitude=forcing.barotropicVelocityAmplitude);
            downsampled = highForcing.forcingWithResolutionOfTransform(sourceTransform);
            testCase.verifyEqual(downsampled.topographicHeight,terrain,AbsTol=1e-12)

            incompatibleTransform = WVTransformBoussinesq([50e3 30e3 2e3],[12 10 7],N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            testCase.verifyError(@()forcing.forcingWithResolutionOfTransform(incompatibleTransform),"WVPseudoTopographicWaveGeneration:IncompatibleDomain")
        end

        function explicitAntialiasConversionPreservesForcing(testCase)
            wvt = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],true);
            terrain = TestWVPseudoTopographicWaveGenerationProduction.topography(wvt);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; -0.02],name="antialias terrain");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            explicitTransform = wvt.waveVortexTransformWithExplicitAntialiasing();
            explicitForcing = explicitTransform.forcingWithName(forcing.name);

            testCase.verifyClass(explicitForcing,"WVPseudoTopographicWaveGeneration")
            testCase.verifyEqual(explicitForcing.topographicHeight,forcing.topographicHeight,AbsTol=1e-12)
            explicitTransform.t = 391;
            diagnostics = TestWVPseudoTopographicWaveGenerationProduction.sourceDiagnostics(explicitTransform,explicitForcing);
            testCase.verifyEqual(diagnostics.F0,zeros(size(explicitTransform.A0)))
            testCase.verifyLessThanOrEqual(diagnostics.qgpvNorm,1e-10)
        end

        function resolutionConversionUsesTargetAdaptiveDamping(testCase)
            source = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],false);
            generation = WVPseudoTopographicWaveGeneration(source,topographicHeight=TestWVPseudoTopographicWaveGenerationProduction.topography(source),barotropicVelocityAmplitude=[0.05; -0.02]);
            target = TestWVPseudoTopographicWaveGenerationProduction.createTransform([12 10 7],true);
            damping = WVAdaptiveDamping(target);
            target.addForcing(damping);
            converted = generation.forcingWithResolutionOfTransform(target);
            target.addForcing(converted);

            [mask,components] = converted.spectralGenerationMask();
            testCase.verifyEqual(components.adaptiveDamping,damping.damp == 0)
            testCase.verifyFalse(any(mask(damping.damp ~= 0),"all"))
            target.t = 391;
            [Fp,Fm] = converted.addSpectralForcing(target,zeros(size(target.Ap)),zeros(size(target.Am)),zeros(size(target.A0)));
            testCase.verifyEqual(Fp(damping.damp ~= 0),zeros(nnz(damping.damp ~= 0),1))
            testCase.verifyEqual(Fm(damping.damp ~= 0),zeros(nnz(damping.damp ~= 0),1))
        end

        function persistenceMetadataAndTransformRoundTrip(testCase)
            requiredProperties = {'topographicHeight','barotropicVelocityAmplitude','frequency','darwinSymbol','rampDuration','startTime','shouldAvoidAdaptiveDamping','maximumForcedHorizontalWavenumber','maximumForcedVerticalMode','name'};
            testCase.verifyEqual(WVPseudoTopographicWaveGeneration.classRequiredPropertyNames(),requiredProperties)
            annotations = WVPseudoTopographicWaveGeneration.classDefinedPropertyAnnotations();
            annotationNames = string({annotations.name});
            testCase.verifyEqual(sort(annotationNames),sort([string(requiredProperties) "barotropicVelocityComponent"]))
            terrainAnnotation = annotations(annotationNames == "topographicHeight");
            velocityAnnotation = annotations(annotationNames == "barotropicVelocityAmplitude");
            testCase.verifyEqual(string(terrainAnnotation.dimensions),["x" "y"])
            testCase.verifyTrue(velocityAnnotation.isComplex)
            testCase.verifyEqual(string(velocityAnnotation.dimensions),"barotropicVelocityComponent")

            wvt = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],false);
            terrain = TestWVPseudoTopographicWaveGenerationProduction.topography(wvt);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.04+0.01i; -0.02+0.015i],frequency=1.31e-4,rampDuration=400,startTime=120,shouldAvoidAdaptiveDamping=false,maximumForcedHorizontalWavenumber=2.5e-4,maximumForcedVerticalMode=3,name="restart terrain");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            wvt.t = 317;
            [originalFp,originalFm,originalF0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));

            path = string(tempname)+".nc";
            cleanup = onCleanup(@()TestWVPseudoTopographicWaveGenerationProduction.deleteFile(path));
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            ncfile.close();
            [restoredTransform,restoredFile] = WVTransform.waveVortexTransformFromFile(path);
            restoredFile.close();
            restoredForcing = restoredTransform.forcingWithName(forcing.name);

            testCase.verifyClass(restoredForcing,"WVPseudoTopographicWaveGeneration")
            testCase.verifyEqual(restoredForcing.topographicHeight,forcing.topographicHeight)
            testCase.verifyEqual(restoredForcing.barotropicVelocityAmplitude,forcing.barotropicVelocityAmplitude)
            testCase.verifyEqual(restoredForcing.frequency,forcing.frequency)
            testCase.verifyEqual(restoredForcing.darwinSymbol,forcing.darwinSymbol)
            testCase.verifyEqual(restoredForcing.rampDuration,forcing.rampDuration)
            testCase.verifyEqual(restoredForcing.startTime,forcing.startTime)
            testCase.verifyEqual(restoredForcing.shouldAvoidAdaptiveDamping,forcing.shouldAvoidAdaptiveDamping)
            testCase.verifyEqual(restoredForcing.maximumForcedHorizontalWavenumber,forcing.maximumForcedHorizontalWavenumber)
            testCase.verifyEqual(restoredForcing.maximumForcedVerticalMode,forcing.maximumForcedVerticalMode)
            testCase.verifyEqual(string(restoredForcing.name),string(forcing.name))
            restoredTransform.t = wvt.t;
            [restoredFp,restoredFm,restoredF0] = restoredForcing.addSpectralForcing(restoredTransform,zeros(size(restoredTransform.Ap)),zeros(size(restoredTransform.Am)),zeros(size(restoredTransform.A0)));
            testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(restoredFp,originalFp),1e-12)
            testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGenerationProduction.relativeError(restoredFm,originalFm),1e-12)
            testCase.verifyEqual(restoredF0,originalF0)
            clear cleanup
        end

        function persistedFrequencyIsAuthoritative(testCase)
            wvt = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],false);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=TestWVPseudoTopographicWaveGenerationProduction.topography(wvt),barotropicVelocityAmplitude=[0.04; -0.01],darwinSymbol="N2",name="authoritative frequency");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);

            path = string(tempname)+".nc";
            cleanup = onCleanup(@()TestWVPseudoTopographicWaveGenerationProduction.deleteFile(path));
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            ncfile.close();
            persistedFrequency = forcing.frequency*(1+1e-6);
            ncwrite(path,"/forcing/frequency",persistedFrequency)

            restoredTransform = WVTransform.waveVortexTransformFromFile(path);
            restoredForcing = restoredTransform.forcingWithName(forcing.name);
            testCase.verifyEqual(restoredForcing.frequency,persistedFrequency)
            testCase.verifyEqual(restoredForcing.darwinSymbol,"N2")

            targetTransform = TestWVPseudoTopographicWaveGenerationProduction.createTransform([12 10 7],false);
            convertedForcing = restoredForcing.forcingWithResolutionOfTransform(targetTransform);
            testCase.verifyEqual(convertedForcing.frequency,persistedFrequency)
            testCase.verifyEqual(convertedForcing.darwinSymbol,"N2")
            clear cleanup
        end

        function adaptiveModelRestartMatchesControl(testCase)
            restartPath = string(tempname)+".nc";
            cleanup = onCleanup(@()TestWVPseudoTopographicWaveGenerationProduction.deleteFile(restartPath));
            [restartModel,checkpointTime,finalTime] = TestWVPseudoTopographicWaveGenerationProduction.modelForRestart();
            [controlModel,~,~] = TestWVPseudoTopographicWaveGenerationProduction.modelForRestart();
            restartModel.createNetCDFFileForModelOutput(restartPath,outputInterval=checkpointTime/2,shouldOverwriteExisting=true);
            restartModel.integrateToTime(checkpointTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            restartModel.closeNetCDFFile();

            controlModel.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel = WVModel.modelFromFile(restartPath);
            resumedModel.setupIntegrator(integratorType="adaptive",absTolerance=1e-12,relTolerance=1e-10);
            resumedModel.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel.closeNetCDFFile();

            testCase.verifyLessThanOrEqual(max(abs(resumedModel.wvt.Ap-controlModel.wvt.Ap),[],"all"),1e-10)
            testCase.verifyLessThanOrEqual(max(abs(resumedModel.wvt.Am-controlModel.wvt.Am),[],"all"),1e-10)
            testCase.verifyEqual(resumedModel.wvt.A0,controlModel.wvt.A0)
            testCase.verifyClass(resumedModel.wvt.forcingWithName("restart model terrain"),"WVPseudoTopographicWaveGeneration")
            clear cleanup
        end

    end

    methods (Static, Access = private)
        function wvt = createTransform(resolution,shouldAntialias)
            wvt = WVTransformBoussinesq([40e3 30e3 2e3],resolution,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=shouldAntialias);
        end

        function terrain = topography(wvt)
            [x,y] = ndgrid(wvt.x,wvt.y);
            terrain = 30*cos(2*pi*x/wvt.Lx)+17*sin(2*pi*y/wvt.Ly)+11*cos(2*pi*(2*x/wvt.Lx+y/wvt.Ly));
        end

        function diagnostics = sourceDiagnostics(wvt,forcing)
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [Fu,Fv,~,Feta] = wvt.transformWaveVortexToUVWEta(Fp,Fm,F0,wvt.t);
            vorticityX = wvt.diffX(Fv);
            vorticityY = wvt.diffY(Fu);
            stretching = wvt.f*wvt.diffZG(Feta);
            qgpvSource = vorticityX-vorticityY-stretching;
            modalPower = 2*sum(wvt.Apm_TE_factor(:).*real(Fp(:).*conj(wvt.Ap(:))+Fm(:).*conj(wvt.Am(:))));
            [~,maskComponents] = forcing.spectralGenerationMask();
            kinematicPressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=wvt.NAp.*wvt.Apt.*maskComponents.effectivePositive+wvt.NAm.*wvt.Amt.*maskComponents.effectiveNegative);
            [~,iBottom] = min(wvt.z);
            bottomPower = mean(kinematicPressure(:,:,iBottom).*forcing.bottomVelocityAtTime(wvt.t),"all");
            diagnostics = struct(F0=F0,modalQGPVSource=wvt.A0_QGPV_factor.*F0,qgpvNorm=norm(qgpvSource(:)),modalPower=modalPower,bottomPower=bottomPower,powerError=abs(modalPower-bottomPower)/max([abs(modalPower) abs(bottomPower) eps]));
        end

        function [model,checkpointTime,finalTime] = modelForRestart()
            wvt = TestWVPseudoTopographicWaveGenerationProduction.createTransform([8 6 5],false);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=TestWVPseudoTopographicWaveGenerationProduction.topography(wvt),barotropicVelocityAmplitude=[0.05; -0.01],rampDuration=100,startTime=0,name="restart model terrain");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            model = WVModel(wvt);
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-12,relTolerance=1e-10);
            checkpointTime = 300;
            finalTime = 600;
        end

        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end

        function deleteFile(path)
            if isfile(path)
                delete(path)
            end
        end
    end
end
