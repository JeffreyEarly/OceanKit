classdef TestForcingLifecycle < matlab.unittest.TestCase
    properties
        tempFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.tempFolder = string(fixture.Folder);
        end
    end

    methods (Test, TestTags="full")
        function registryMutationsAreIdentityBasedAndAtomic(testCase)
            wvt = TestForcingLifecycle.constantTransform([8 6 5]);
            wvt.removeAllForcing();
            spatialLate = WVTestForcing(wvt,"spatial late",WVForcingType("HydrostaticSpatial"),uint8(200),1);
            spatialFirst = WVTestForcing(wvt,"spatial first",WVForcingType("HydrostaticSpatial"),uint8(10),2);
            spatialTie = WVTestForcing(wvt,"spatial tie",WVForcingType("HydrostaticSpatial"),uint8(200),3);
            spectral = WVTestForcing(wvt,"spectral",WVForcingType("Spectral"),uint8(1),4);
            amplitude = WVTestForcing(wvt,"amplitude",WVForcingType("SpectralAmplitude"),uint8(1),5);

            wvt.addForcing([spatialLate spatialFirst spatialTie spectral amplitude]);
            testCase.verifyEqual(wvt.forcingNames,["spatial first"; "spatial late"; "spatial tie"; "spectral"; "amplitude"])
            expected = wvt.forcing;

            wvt.addForcing(spatialLate);
            TestForcingLifecycle.verifyIdentity(testCase,wvt.forcing,expected)

            replacement = WVTestForcing(wvt,"spatial late",WVForcingType("HydrostaticSpatial"),uint8(200),6);
            wvt.addForcing(replacement);
            testCase.verifyTrue(wvt.forcingWithName("spatial late") == replacement)
            testCase.verifyEqual(wvt.forcingNames,["spatial first"; "spatial tie"; "spatial late"; "spectral"; "amplitude"])

            firstDuplicate = WVTestForcing(wvt,"batch duplicate",WVForcingType("HydrostaticSpatial"),uint8(200),7);
            lastDuplicate = WVTestForcing(wvt,"batch duplicate",WVForcingType("HydrostaticSpatial"),uint8(200),8);
            wvt.addForcing([firstDuplicate lastDuplicate]);
            testCase.verifyTrue(wvt.forcingWithName("batch duplicate") == lastDuplicate)

            beforeFailure = wvt.forcing;
            valid = WVTestForcing(wvt,"would be valid",WVForcingType("Spectral"),uint8(20),9);
            unsupported = WVTestForcing(wvt,"unsupported",WVForcingType("PVSpectral"),uint8(20),10);
            testCase.verifyError(@()wvt.addForcing([valid unsupported]),'')
            TestForcingLifecycle.verifyIdentity(testCase,wvt.forcing,beforeFailure)
            testCase.verifyFalse(wvt.hasForcingWithName("would be valid"))

            otherTransform = TestForcingLifecycle.constantTransform([8 6 5]);
            foreign = WVTestForcing(otherTransform,"spatial first",WVForcingType("HydrostaticSpatial"),uint8(10),11);
            testCase.verifyError(@()wvt.addForcing(foreign),'')
            testCase.verifyError(@()wvt.removeForcing(foreign),'')
            testCase.verifyError(@()wvt.forcingWithName("missing"),'')
            TestForcingLifecycle.verifyIdentity(testCase,wvt.forcing,beforeFailure)

            testCase.verifyError(@()wvt.setForcing([valid unsupported]),'')
            TestForcingLifecycle.verifyIdentity(testCase,wvt.forcing,beforeFailure)

            wvt.removeForcing([replacement replacement]);
            testCase.verifyFalse(wvt.hasForcingWithName("spatial late"))
            testCase.verifyError(@()wvt.removeForcing(replacement),'')

            columnForcing = reshape([WVTestForcing(wvt,"column one",WVForcingType("Spectral"),uint8(30),1) ...
                WVTestForcing(wvt,"column two",WVForcingType("Spectral"),uint8(40),2)],[],1);
            wvt.addForcing(columnForcing);
            testCase.verifyTrue(all(wvt.hasForcingWithName("column one","column two")))
        end

        function suppliedForcingCoversCompatibleLifecycleAndFluxes(testCase)
            seedRandomNumberGenerator(testCase,3811);
            wvt = TestForcingLifecycle.boussinesqTransform([8 6 5]);
            wvt.removeAllForcing();
            forces = TestForcingLifecycle.waveForcingInventory(wvt);
            wvt.addForcing(forces);

            expectedClasses = ["WVNonlinearAdvection" "WVAdaptiveDamping" "WVAntialiasing" ...
                "WVHorizontalDamping" "WVVerticalDamping" "WVVerticalDiffusivity" ...
                "WVFixedAmplitudeForcing" "WVBottomFrictionLinear" ...
                "WVBottomFrictionQuadratic" "WVBetaPlanePVAdvection" ...
                "WVPseudoTopographicWaveGeneration"];
            testCase.verifyEqual(sort(string(arrayfun(@class,wvt.forcing,UniformOutput=false))),sort(expectedClasses))
            for iForce = 1:length(forces)
                testCase.verifyTrue(wvt.forcingWithName(forces(iForce).name) == forces(iForce))
            end

            wvt.initWithRandomFlow(uvMax=0.01);
            [Fp,Fm,F0] = wvt.nonlinearFlux;
            TestForcingLifecycle.verifyFinite(testCase,Fp,Fm,F0)

            adaptive = wvt.forcingWithName("adaptive damping");
            [adaptiveFp,adaptiveFm,adaptiveF0] = adaptive.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(adaptiveFp,wvt.uvMax*adaptive.damp.*wvt.Ap)
            testCase.verifyEqual(adaptiveFm,wvt.uvMax*adaptive.damp.*wvt.Am)
            testCase.verifyEqual(adaptiveF0,wvt.uvMax*adaptive.damp.*wvt.A0)

            antialias = wvt.forcingWithName("antialias filter");
            onesAp = ones(size(wvt.Ap));
            [filteredAp,filteredAm,filteredA0] = antialias.addSpectralForcing(wvt,onesAp,onesAp,onesAp);
            testCase.verifyEqual(filteredAp,double(~antialias.M))
            testCase.verifyEqual(filteredAm,double(~antialias.M))
            testCase.verifyEqual(filteredA0,double(~antialias.M))

            fixed = wvt.forcingWithName("fixed modes");
            [fixedAp,fixedAm,fixedA0] = fixed.setSpectralAmplitude(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(fixedAp(fixed.Ap_indices),fixed.Apbar)
            testCase.verifyEqual(fixedAm(fixed.Am_indices),fixed.Ambar)
            testCase.verifyEqual(fixedA0(fixed.A0_indices),fixed.A0bar)

            vertical = wvt.forcingWithName("vertical diffusivity");
            [~,~,~,Feta] = vertical.addNonhydrostaticSpatialForcing(wvt,zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize),zeros(wvt.spatialMatrixSize));
            testCase.verifyEqual(Feta,vertical.kappa_z*wvt.diffZG(wvt.eta,n=2),AbsTol=1e-12)

            TestForcingLifecycle.verifyBottomFriction(testCase,wvt)
            TestForcingLifecycle.verifyNonlinearAdvection(testCase,wvt)

            qg = TestForcingLifecycle.barotropicTransform([8 6]);
            qg.removeAllForcing();
            qgForces = [WVNonlinearAdvection(qg) WVBottomFrictionLinear(qg,r=2e-7) WVBottomFrictionQuadratic(qg,Cd=1.5e-3) ...
                WVBetaPlanePVAdvection(qg) WVAdaptiveDamping(qg) WVAntialiasing(qg,Nj=1) ...
                WVFixedAmplitudeForcing(qg,name="fixed qg",A0_indices=uint64(2),A0bar=3e-3)];
            qg.addForcing(qgForces);
            barotropicDiffusivity = WVVerticalDiffusivity(qg,kappa_z=7e-6);
            testCase.verifyError(@()qg.addForcing(barotropicDiffusivity),'')
            qg.A0 = complex(randn(size(qg.A0)),randn(size(qg.A0))).*qg.geostrophicComponent.maskA0;
            F0qg = qg.nonlinearFlux;
            TestForcingLifecycle.verifyFinite(testCase,F0qg)
            beta = qg.forcingWithName("beta-plane advection of qgpv");
            testCase.verifyEqual(beta.addPotentialVorticitySpatialForcing(qg,zeros(qg.spatialMatrixSize)),-qg.beta*qg.v)

            stratifiedQG = WVTransformStratifiedQG([40e3 30e3 2e3],[8 6 5],N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            stratifiedQG.removeAllForcing();
            qgDiffusivity = WVVerticalDiffusivity(stratifiedQG,kappa_z=7e-6,shouldForceMeanDensityAnomaly=false);
            stratifiedQG.addForcing(qgDiffusivity);
            stratifiedQG.A0 = complex(randn(size(stratifiedQG.A0)),randn(size(stratifiedQG.A0))).*stratifiedQG.geostrophicComponent.maskA0;
            TestForcingLifecycle.verifyFinite(testCase,stratifiedQG.nonlinearFlux)
        end

        function resolutionConversionPreservesConfigurationAndSelections(testCase)
            source = TestForcingLifecycle.boussinesqTransform([8 6 5]);
            source.removeAllForcing();
            vertical = WVVerticalDiffusivity(source,kappa_z=7.25e-6,shouldForceMeanDensityAnomaly=false);
            antialias = WVAntialiasing(source,Nj=3);
            indices = uint64([2; 4; 6]);
            fixed = WVFixedAmplitudeForcing(source,name="selected modes",Ap_indices=indices,Apbar=[0; 2+3i; -4i],Am_indices=indices,Ambar=[1-2i; 0; 5],A0_indices=indices,A0bar=[0; -3; 2i]);
            source.setForcing([vertical antialias fixed]);

            target = source.waveVortexTransformWithResolution([12 10 7]);
            convertedVertical = target.forcingWithName("vertical diffusivity");
            testCase.verifyEqual(convertedVertical.kappa_z,vertical.kappa_z)
            testCase.verifyEqual(convertedVertical.shouldForceMeanDensityAnomaly,vertical.shouldForceMeanDensityAnomaly)
            testCase.verifyEqual(target.forcingWithName("antialias filter").Nj,3)
            TestForcingLifecycle.verifyConvertedFixedAmplitude(testCase,source,target,fixed)

            downsampled = target.waveVortexTransformWithResolution([6 4 3]);
            testCase.verifyEqual(downsampled.forcingWithName("vertical diffusivity").kappa_z,vertical.kappa_z)
            testCase.verifyEqual(downsampled.forcingWithName("antialias filter").Nj,3)
            TestForcingLifecycle.verifyConvertedFixedAmplitude(testCase,target,downsampled,target.forcingWithName("selected modes"))

            inventorySource = TestForcingLifecycle.boussinesqTransform([8 6 5]);
            inventorySource.removeAllForcing();
            inventorySource.addForcing(TestForcingLifecycle.waveForcingInventory(inventorySource));
            convertedInventory = inventorySource.waveVortexTransformWithResolution([12 10 7]);
            testCase.verifyEqual(convertedInventory.forcingNames,inventorySource.forcingNames)
            testCase.verifyEqual(string(arrayfun(@class,convertedInventory.forcing,UniformOutput=false)),string(arrayfun(@class,inventorySource.forcing,UniformOutput=false)))
            testCase.verifyEqual(convertedInventory.forcingWithName("horizontal scalar diffusivity").nu,0.125)
            testCase.verifyEqual(convertedInventory.forcingWithName("vertical scalar diffusivity").kappa,7.5e-5)
            testCase.verifyEqual(convertedInventory.forcingWithName("linear bottom friction").r,2e-7)
            testCase.verifyEqual(convertedInventory.forcingWithName("quadratic bottom friction").Cd,1.5e-3)
            testCase.verifyEqual(convertedInventory.forcingWithName("antialias filter").Nj,3)
            convertedInventory.initWithRandomFlow(uvMax=0.01);
            [convertedFp,convertedFm,convertedF0] = convertedInventory.nonlinearFlux;
            TestForcingLifecycle.verifyFinite(testCase,convertedFp,convertedFm,convertedF0)
        end

        function composedForcingMatchesIndependentContributions(testCase)
            seedRandomNumberGenerator(testCase,3812);
            wvt = TestForcingLifecycle.constantTransform([8 6 5]);
            wvt.removeAllForcing();
            horizontal = WVHorizontalDamping(wvt,nu=0.125,kappa=2.5e-5);
            vertical = WVVerticalDamping(wvt,nu=0.375,kappa=7.5e-5);
            adaptive = WVAdaptiveDamping(wvt);
            [x,y] = ndgrid(wvt.x,wvt.y);
            generation = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=20*cos(2*pi*x/wvt.Lx)+10*sin(2*pi*y/wvt.Ly),barotropicVelocityAmplitude=[0.04;-0.02],frequency=1.31e-4);
            wvt.addForcing([horizontal vertical adaptive generation]);
            wvt.initWithRandomFlow(uvMax=0.01);
            wvt.t = 217;

            zeroField = zeros(wvt.spatialMatrixSize);
            [horizontalFu,horizontalFv,horizontalFw,horizontalFeta] = horizontal.addNonhydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField,zeroField);
            [verticalFu,verticalFv,verticalFw,verticalFeta] = vertical.addNonhydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField,zeroField);
            [expectedFp,expectedFm,expectedF0] = wvt.transformUVWEtaToWaveVortex(horizontalFu+verticalFu,horizontalFv+verticalFv,horizontalFw+verticalFw,horizontalFeta+verticalFeta);
            [expectedFp,expectedFm,expectedF0] = adaptive.addSpectralForcing(wvt,expectedFp,expectedFm,expectedF0);
            [expectedFp,expectedFm,expectedF0] = generation.addSpectralForcing(wvt,expectedFp,expectedFm,expectedF0);
            [actualFp,actualFm,actualF0] = wvt.nonlinearFlux;
            testCase.verifyEqual(actualFp,expectedFp,AbsTol=1e-12)
            testCase.verifyEqual(actualFm,expectedFm,AbsTol=1e-12)
            testCase.verifyEqual(actualF0,expectedF0,AbsTol=1e-12)

            firstFixed = WVFixedAmplitudeForcing(wvt,name="first fixed",Ap_indices=uint64(2),Apbar=1+2i);
            secondFixed = WVFixedAmplitudeForcing(wvt,name="second fixed",Am_indices=uint64(4),Ambar=-3i,A0_indices=uint64(6),A0bar=2);
            wvt.addForcing([firstFixed secondFixed]);
            wvt.restoreForcingAmplitudes();
            testCase.verifyEqual(wvt.Ap(2),1+2i)
            testCase.verifyEqual(wvt.Am(4),-3i)
            testCase.verifyEqual(wvt.A0(6),2)

            qg = TestForcingLifecycle.barotropicTransform([8 6]);
            qg.removeAllForcing();
            nonlinear = WVNonlinearAdvection(qg);
            linearBottom = WVBottomFrictionLinear(qg,r=2e-7);
            beta = WVBetaPlanePVAdvection(qg);
            qg.addForcing([nonlinear linearBottom beta]);
            qg.A0 = complex(randn(size(qg.A0)),randn(size(qg.A0))).*qg.geostrophicComponent.maskA0;
            zeroPV = zeros(qg.spatialMatrixSize);
            independentPV = nonlinear.addPotentialVorticitySpatialForcing(qg,zeroPV) ...
                +linearBottom.addPotentialVorticitySpatialForcing(qg,zeroPV) ...
                +beta.addPotentialVorticitySpatialForcing(qg,zeroPV);
            testCase.verifyEqual(qg.nonlinearFlux,qg.transformQGPVToWaveVortex(independentPV),AbsTol=1e-12)
        end

        function netCDFRoundTripPreservesForcingRegistryAndBehavior(testCase)
            wvt = TestForcingLifecycle.boussinesqTransform([8 6 5]);
            wvt.removeAllForcing();
            wvt.addForcing(TestForcingLifecycle.waveForcingInventory(wvt));
            wvt.initWithRandomFlow(uvMax=0.01);
            wvt.t = 317;
            [Fp,Fm,F0] = wvt.nonlinearFlux;
            expectedNames = wvt.forcingNames;

            path = fullfile(testCase.tempFolder,"stable-forcing-round-trip.nc");
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            cleanup = onCleanup(@()TestForcingLifecycle.closeIfOpen(ncfile));
            ncfile.close();
            clear cleanup

            [restored,restoredFile] = WVTransform.waveVortexTransformFromFile(path);
            restoredCleanup = onCleanup(@()TestForcingLifecycle.closeIfOpen(restoredFile));
            testCase.verifyEqual(restored.forcingNames,expectedNames)
            testCase.verifyEqual(restored.forcingWithName("antialias filter").Nj,3)
            testCase.verifyEqual(restored.forcingWithName("vertical diffusivity").kappa_z,7e-6)
            testCase.verifyFalse(restored.forcingWithName("vertical diffusivity").shouldForceMeanDensityAnomaly)
            restored.t = wvt.t;
            [restoredFp,restoredFm,restoredF0] = restored.nonlinearFlux;
            testCase.verifyEqual(restoredFp,Fp,AbsTol=1e-12)
            testCase.verifyEqual(restoredFm,Fm,AbsTol=1e-12)
            testCase.verifyEqual(restoredF0,F0,AbsTol=1e-12)
            restoredFile.close();
            clear restoredCleanup

            reopened = NetCDFFile(path,shouldReadOnly=true);
            reopened.close();

            qg = TestForcingLifecycle.barotropicTransform([8 6]);
            qg.removeAllForcing();
            qg.addForcing([WVNonlinearAdvection(qg) WVBottomFrictionLinear(qg,r=2e-7) ...
                WVBottomFrictionQuadratic(qg,Cd=1.5e-3) WVBetaPlanePVAdvection(qg) ...
                WVAdaptiveDamping(qg) WVAntialiasing(qg,Nj=1) ...
                WVFixedAmplitudeForcing(qg,name="fixed qg",A0_indices=uint64(2),A0bar=3e-3)]);
            qg.A0 = complex(reshape(1:numel(qg.A0),size(qg.A0)),reshape(numel(qg.A0):-1:1,size(qg.A0))).*qg.geostrophicComponent.maskA0;
            expectedQGFlux = qg.nonlinearFlux;
            qgPath = fullfile(testCase.tempFolder,"qg-forcing-round-trip.nc");
            qgFile = qg.writeToFile(qgPath,shouldOverwriteExisting=true);
            qgCleanup = onCleanup(@()TestForcingLifecycle.closeIfOpen(qgFile));
            qgFile.close();
            clear qgCleanup
            [restoredQG,restoredQGFile] = WVTransform.waveVortexTransformFromFile(qgPath);
            restoredQGCleanup = onCleanup(@()TestForcingLifecycle.closeIfOpen(restoredQGFile));
            testCase.verifyEqual(restoredQG.forcingNames,qg.forcingNames)
            testCase.verifyEqual(restoredQG.forcingWithName("antialias filter").Nj,1)
            testCase.verifyEqual(restoredQG.forcingWithName("linear bottom friction").r,2e-7)
            testCase.verifyEqual(restoredQG.forcingWithName("quadratic bottom friction").Cd,1.5e-3)
            testCase.verifyEqual(restoredQG.nonlinearFlux,expectedQGFlux,AbsTol=1e-12)
            restoredQGFile.close();
            clear restoredQGCleanup
        end
    end

    methods (Static, Access=private)
        function wvt = constantTransform(resolution)
            wvt = WVTransformConstantStratification([40e3 30e3 2e3],resolution,N0=5.2e-3,latitude=45,isHydrostatic=false,shouldAntialias=false);
        end

        function wvt = boussinesqTransform(resolution)
            wvt = WVTransformBoussinesq([40e3 30e3 2e3],resolution,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
        end

        function wvt = barotropicTransform(resolution)
            wvt = WVTransformBarotropicQG([40e3 30e3],resolution,latitude=45,shouldAntialias=false);
        end

        function forcing = waveForcingInventory(wvt)
            [x,y] = ndgrid(wvt.x,wvt.y);
            topography = 20*cos(2*pi*x/wvt.Lx)+10*sin(2*pi*y/wvt.Ly);
            forcing = [WVNonlinearAdvection(wvt) WVAdaptiveDamping(wvt) WVAntialiasing(wvt,Nj=3) ...
                WVHorizontalDamping(wvt,nu=0.125,kappa=2.5e-5) WVVerticalDamping(wvt,nu=0.375,kappa=7.5e-5) ...
                WVVerticalDiffusivity(wvt,kappa_z=7e-6,shouldForceMeanDensityAnomaly=false) ...
                WVFixedAmplitudeForcing(wvt,name="fixed modes",Ap_indices=uint64([2;4]),Apbar=[0;2+3i],Am_indices=uint64([2;4]),Ambar=[1-2i;0],A0_indices=uint64([2;4]),A0bar=[0;-3]) ...
                WVBottomFrictionLinear(wvt,r=2e-7) WVBottomFrictionQuadratic(wvt,Cd=1.5e-3) ...
                WVBetaPlanePVAdvection(wvt) WVPseudoTopographicWaveGeneration(wvt,topographicHeight=topography,barotropicVelocityAmplitude=[0.04;-0.02],frequency=1.31e-4,name="terrain generation")];
        end

        function verifyBottomFriction(testCase,wvt)
            linear = wvt.forcingWithName("linear bottom friction");
            quadratic = wvt.forcingWithName("quadratic bottom friction");
            zerosField = zeros(wvt.spatialMatrixSize);
            [linearFu,linearFv] = linear.addNonhydrostaticSpatialForcing(wvt,zerosField,zerosField,zerosField,zerosField);
            testCase.verifyEqual(linearFu(:,:,1),-linear.r_scaled*wvt.u(:,:,1),AbsTol=1e-24)
            testCase.verifyEqual(linearFv(:,:,1),-linear.r_scaled*wvt.v(:,:,1),AbsTol=1e-24)
            testCase.verifyEqual(linearFu(:,:,2:end),zeros(size(linearFu(:,:,2:end))))
            [quadraticFu,quadraticFv] = quadratic.addNonhydrostaticSpatialForcing(wvt,zerosField,zerosField,zerosField,zerosField);
            speed = hypot(wvt.u(:,:,1),wvt.v(:,:,1));
            testCase.verifyEqual(quadraticFu(:,:,1),-quadratic.cd*wvt.u(:,:,1).*speed,AbsTol=1e-24)
            testCase.verifyEqual(quadraticFv(:,:,1),-quadratic.cd*wvt.v(:,:,1).*speed,AbsTol=1e-24)
        end

        function verifyNonlinearAdvection(testCase,wvt)
            nonlinear = wvt.forcingWithName("nonlinear advection");
            zerosField = zeros(wvt.spatialMatrixSize);
            [Fu,Fv,Fw,Feta] = nonlinear.addNonhydrostaticSpatialForcing(wvt,zerosField,zerosField,zerosField,zerosField);
            expectedFu = -(wvt.u.*wvt.diffX(wvt.u)+wvt.v.*wvt.diffY(wvt.u)+wvt.w.*wvt.diffZF(wvt.u));
            expectedFv = -(wvt.u.*wvt.diffX(wvt.v)+wvt.v.*wvt.diffY(wvt.v)+wvt.w.*wvt.diffZF(wvt.v));
            expectedFw = -(wvt.u.*wvt.diffX(wvt.w)+wvt.v.*wvt.diffY(wvt.w)+wvt.w.*wvt.diffZG(wvt.w));
            expectedFeta = -(wvt.u.*wvt.diffX(wvt.eta)+wvt.v.*wvt.diffY(wvt.eta)+wvt.w.*(wvt.diffZG(wvt.eta)+wvt.eta.*nonlinear.dLnN2));
            testCase.verifyEqual(Fu,expectedFu)
            testCase.verifyEqual(Fv,expectedFv)
            testCase.verifyEqual(Fw,expectedFw)
            testCase.verifyEqual(Feta,expectedFeta)
        end

        function verifyConvertedFixedAmplitude(testCase,source,target,sourceForcing)
            targetForcing = target.forcingWithName(sourceForcing.name);
            coefficientNames = ["Ap" "Am" "A0"];
            for coefficientName = coefficientNames
                sourceValues = complex(zeros(source.spectralMatrixSize));
                sourceIndices = sourceForcing.(coefficientName+"_indices");
                sourceValues(sourceIndices) = sourceForcing.(coefficientName+"bar");
                expectedValues = source.spectralVariableWithResolution(target,sourceValues);
                sourceMask = zeros(source.spectralMatrixSize);
                sourceMask(sourceIndices) = 1;
                expectedMask = source.spectralVariableWithResolution(target,sourceMask);
                expectedIndices = find(expectedMask);
                testCase.verifyEqual(targetForcing.(coefficientName+"_indices"),uint64(expectedIndices))
                testCase.verifyEqual(targetForcing.(coefficientName+"bar"),expectedValues(expectedIndices))
            end
        end

        function verifyIdentity(testCase,actual,expected)
            testCase.verifyEqual(length(actual),length(expected))
            testCase.verifyTrue(all(actual == expected))
        end

        function verifyFinite(testCase,varargin)
            for iValue = 1:length(varargin)
                testCase.verifyTrue(all(isfinite(varargin{iValue}),"all"))
            end
        end


        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
