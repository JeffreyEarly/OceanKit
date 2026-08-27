classdef TestNarrowBandGeostrophicForcing < matlab.unittest.TestCase
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
        function seededSubclassMatchesCompatibilityHelper(testCase)
            seedRandomNumberGenerator(testCase,4201);
            rng(4201,"twister")
            wvt = TestNarrowBandGeostrophicForcing.stratifiedTransform([8 8 5]);
            k_r = wvt.dk;
            k_f = 2*wvt.dk;
            force = WVNarrowBandGeostrophicForcing(wvt,name="direct forcing",k_r=k_r,k_f=k_f,j_f=1,u_rms=0.04,initialPV="full-spectrum");
            A0 = wvt.A0;
            nextRandomValues = rand(1,5);

            rng(4201,"twister")
            legacyWvt = TestNarrowBandGeostrophicForcing.stratifiedTransform([8 8 5]);
            legacy = WVFixedAmplitudeForcing(legacyWvt,name="direct forcing");
            [legacySpectrum,legacyR] = legacy.setNarrowBandGeostrophicForcing(k_r=k_r,k_f=k_f,j_f=1,u_rms=0.04,initialPV="full-spectrum");
            legacyNextRandomValues = rand(1,5);

            testCase.verifyEqual(legacyWvt.A0,A0)
            testCase.verifyEqual(legacy.A0_indices,force.A0_indices)
            testCase.verifyEqual(legacy.A0bar,force.A0bar)
            testCase.verifyEqual(legacyR,force.r)
            testCase.verifyEqual(legacyNextRandomValues,nextRandomValues)
            spectrumWavenumbers = [0.5*k_r k_r 1.5*k_r k_f 2*k_f];
            testCase.verifyEqual(legacySpectrum(spectrumWavenumbers),force.modelSpectrum(spectrumWavenumbers))

            F = wvt.FinvMatrix;
            expectedR = 2.25*abs(F(end,2)/F(1,2))*force.u_rms*k_r;
            testCase.verifyEqual(force.r,expectedR)
            testCase.verifyEqual(force.k_r,k_r)
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            expectedIndices = find(wvt.Kh > k_f-deltaK/2 & wvt.Kh < k_f+deltaK/2 & wvt.J == 1);
            testCase.verifyEqual(force.A0_indices,uint64(expectedIndices))
            testCase.verifyEqual(force.A0bar,wvt.A0(expectedIndices))

            rng(4202,"twister")
            barotropic = TestNarrowBandGeostrophicForcing.barotropicTransform([8 8]);
            u_rms = 0.05;
            specifiedR = 0.0225*u_rms*barotropic.dk;
            barotropicForce = WVNarrowBandGeostrophicForcing(barotropic,name="barotropic forcing",r=specifiedR,k_f=2*barotropic.dk,u_rms=u_rms,initialPV="narrow-band");
            barotropicA0 = barotropic.A0;
            barotropicNextRandomValues = rand(1,5);

            rng(4202,"twister")
            legacyBarotropic = TestNarrowBandGeostrophicForcing.barotropicTransform([8 8]);
            legacyBarotropicForce = WVFixedAmplitudeForcing(legacyBarotropic,name="barotropic forcing");
            [legacyBarotropicSpectrum,legacyBarotropicR] = legacyBarotropicForce.setNarrowBandGeostrophicForcing(r=specifiedR,k_f=2*legacyBarotropic.dk,u_rms=u_rms,initialPV="narrow-band");
            legacyBarotropicNextRandomValues = rand(1,5);

            testCase.verifyEqual(legacyBarotropic.A0,barotropicA0)
            testCase.verifyEqual(legacyBarotropicForce.A0_indices,barotropicForce.A0_indices)
            testCase.verifyEqual(legacyBarotropicForce.A0bar,barotropicForce.A0bar)
            testCase.verifyEqual(legacyBarotropicR,barotropicForce.r)
            testCase.verifyEqual(barotropicForce.r,specifiedR)
            testCase.verifyEqual(barotropicForce.k_r,specifiedR/(0.0225*u_rms))
            testCase.verifyEqual(legacyBarotropicSpectrum(spectrumWavenumbers),barotropicForce.modelSpectrum(spectrumWavenumbers))
            testCase.verifyEqual(legacyBarotropicNextRandomValues,barotropicNextRandomValues)
        end

        function netCDFRoundTripUsesCanonicalRestoration(testCase)
            seedRandomNumberGenerator(testCase,4203);
            wvt = TestNarrowBandGeostrophicForcing.stratifiedTransform([8 8 5]);
            wvt.removeAllForcing();
            force = WVNarrowBandGeostrophicForcing(wvt,name="restart forcing",k_r=wvt.dk,k_f=2*wvt.dk,j_f=1,u_rms=0.03,initialPV="narrow-band");
            wvt.addForcing(force);
            expectedA0 = wvt.A0;
            forcingInput = complex(reshape(1:numel(wvt.A0),size(wvt.A0)),reshape(numel(wvt.A0):-1:1,size(wvt.A0)));
            expectedTendency = force.setPotentialVorticitySpectralForcing(wvt,forcingInput);
            expectedAmplitude = force.setPotentialVorticitySpectralAmplitude(wvt,zeros(size(wvt.A0)));

            filePath = fullfile(testCase.tempFolder,"narrow-band-forcing.nc");
            ncfile = wvt.writeToFile(filePath,shouldOverwriteExisting=true);
            cleanup = onCleanup(@()TestNarrowBandGeostrophicForcing.closeIfOpen(ncfile));
            ncfile.close();
            clear cleanup

            rng(4204,"twister")
            randomStateBeforeRestore = rng;
            [restored,restoredFile] = WVTransform.waveVortexTransformFromFile(filePath);
            restoredCleanup = onCleanup(@()TestNarrowBandGeostrophicForcing.closeIfOpen(restoredFile));
            randomStateAfterRestore = rng;
            restoredForce = restored.forcingWithName("restart forcing");

            testCase.verifyClass(restoredForce,"WVNarrowBandGeostrophicForcing")
            TestNarrowBandGeostrophicForcing.verifyConfiguration(testCase,restoredForce,force)
            testCase.verifyEqual(restoredForce.A0_indices,force.A0_indices)
            testCase.verifyEqual(restoredForce.A0bar,force.A0bar)
            testCase.verifyEqual(restored.A0,expectedA0)
            testCase.verifyEqual(randomStateAfterRestore,randomStateBeforeRestore)
            testCase.verifyEqual(restoredForce.setPotentialVorticitySpectralForcing(restored,forcingInput),expectedTendency)
            testCase.verifyEqual(restoredForce.setPotentialVorticitySpectralAmplitude(restored,zeros(size(restored.A0))),expectedAmplitude)
            testCase.verifyEqual(restoredForce.modelSpectrum([0.5*force.k_r force.k_r force.k_f 2*force.k_f]),force.modelSpectrum([0.5*force.k_r force.k_r force.k_f 2*force.k_f]))
            restoredFile.close();
            clear restoredCleanup
        end

        function resolutionConversionPreservesStateWithoutInitialization(testCase)
            seedRandomNumberGenerator(testCase,4205);
            source = TestNarrowBandGeostrophicForcing.stratifiedTransform([8 8 5]);
            source.removeAllForcing();
            source.A0 = complex(reshape(1:numel(source.A0),size(source.A0)),reshape(numel(source.A0):-1:1,size(source.A0))).*source.geostrophicComponent.maskA0;
            force = WVNarrowBandGeostrophicForcing(source,name="converted forcing",k_r=source.dk,k_f=2*source.dk,j_f=1,u_rms=0.025,initialPV="none");

            target = TestNarrowBandGeostrophicForcing.stratifiedTransform([12 10 7]);
            target.removeAllForcing();
            targetA0 = complex(reshape(1:numel(target.A0),size(target.A0)),reshape(numel(target.A0):-1:1,size(target.A0))).*target.geostrophicComponent.maskA0;
            target.A0 = targetA0;
            randomStateBeforeConversion = rng;
            converted = force.forcingWithResolutionOfTransform(target);
            randomStateAfterConversion = rng;

            sourceSelection = zeros(source.spectralMatrixSize);
            sourceSelection(force.A0_indices) = 1;
            targetSelection = source.spectralVariableWithResolution(target,sourceSelection);
            expectedIndices = find(targetSelection);
            sourceValues = complex(zeros(source.spectralMatrixSize));
            sourceValues(force.A0_indices) = force.A0bar;
            targetValues = source.spectralVariableWithResolution(target,sourceValues);

            testCase.verifyClass(converted,"WVNarrowBandGeostrophicForcing")
            TestNarrowBandGeostrophicForcing.verifyConfiguration(testCase,converted,force)
            testCase.verifyEqual(converted.A0_indices,uint64(expectedIndices))
            testCase.verifyEqual(converted.A0bar,targetValues(expectedIndices))
            testCase.verifyEqual(target.A0,targetA0)
            testCase.verifyEqual(randomStateAfterConversion,randomStateBeforeConversion)
        end
    end

    methods (Static, Access = private)
        function wvt = stratifiedTransform(resolution)
            wvt = WVTransformStratifiedQG([40e3 40e3 2e3],resolution,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
        end

        function wvt = barotropicTransform(resolution)
            wvt = WVTransformBarotropicQG([40e3 40e3],resolution,latitude=45,shouldAntialias=false);
        end

        function verifyConfiguration(testCase,actual,expected)
            testCase.verifyEqual(actual.name,expected.name)
            testCase.verifyEqual(actual.r,expected.r)
            testCase.verifyEqual(actual.k_r,expected.k_r)
            testCase.verifyEqual(actual.k_f,expected.k_f)
            testCase.verifyEqual(actual.j_f,expected.j_f)
            testCase.verifyEqual(actual.u_rms,expected.u_rms)
            testCase.verifyEqual(actual.initialPV,expected.initialPV)
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
