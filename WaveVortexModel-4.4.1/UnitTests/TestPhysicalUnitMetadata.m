classdef TestPhysicalUnitMetadata < matlab.unittest.TestCase
    properties
        repositoryRoot
        temporaryFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function stableTransformsUseCanonicalUnits(testCase)
            transforms = TestPhysicalUnitMetadata.stableTransforms();
            for iTransform = 1:numel(transforms)
                testCase.verifyCanonicalAnnotations(transforms{iTransform}.propertyAnnotations(),class(transforms{iTransform}));
            end

            boussinesq = transforms{4};
            testCase.verifyAnnotationUnit(boussinesq.propertyAnnotations(),"K2unique","rad2 m-2");
            testCase.verifyAnnotationUnit(boussinesq.propertyAnnotations(),"Ap","m s-1");
            testCase.verifyAnnotationUnit(boussinesq.propertyAnnotations(),"A0","m2 s-1");
            operationSource = string(fileread(fullfile(testCase.repositoryRoot,"@WVTransform","defaultOperations.m")));
            testCase.verifyMatches(operationSource,"WVVariableAnnotation\('F0',[^\n]*'m2 s-2'");
        end

        function suppliedSubsystemsUseCanonicalUnits(testCase)
            coefficientAnnotations = WVCoefficients.classDefinedPropertyAnnotations();
            tracerAnnotations = WVTracer.classDefinedPropertyAnnotations();
            annotations = [ ...
                WVAdaptiveDamping.classDefinedPropertyAnnotations(), ...
                WVAntialiasing.classDefinedPropertyAnnotations(), ...
                WVBetaPlanePVAdvection.classDefinedPropertyAnnotations(), ...
                WVBottomFrictionLinear.classDefinedPropertyAnnotations(), ...
                WVBottomFrictionQuadratic.classDefinedPropertyAnnotations(), ...
                WVFixedAmplitudeForcing.classDefinedPropertyAnnotations(), ...
                WVHorizontalDamping.classDefinedPropertyAnnotations(), ...
                WVNonlinearAdvection.classDefinedPropertyAnnotations(), ...
                WVPseudoTopographicWaveGeneration.classDefinedPropertyAnnotations(), ...
                WVThermalDamping.classDefinedPropertyAnnotations(), ...
                WVVerticalDamping.classDefinedPropertyAnnotations(), ...
                WVVerticalDiffusivity.classDefinedPropertyAnnotations(), ...
                coefficientAnnotations, ...
                WVEulerianFields.classDefinedPropertyAnnotations(), ...
                WVLagrangianParticles.classDefinedPropertyAnnotations(), ...
                WVMooring.classDefinedPropertyAnnotations(), ...
                tracerAnnotations, ...
                WVModelOutputGroupEvenlySpaced.classDefinedPropertyAnnotations()];
            testCase.verifyCanonicalAnnotations(annotations,"supplied subsystems");
            testCase.verifyAnnotationUnit(coefficientAnnotations,"absTolerance","m2 s-1");
            testCase.verifyAnnotationUnit(tracerAnnotations,"absTolerance","");
            testCase.verifyAnnotationUnit(annotations,"Apbar","m s-1");
            testCase.verifyAnnotationUnit(annotations,"Ambar","m s-1");
            testCase.verifyAnnotationUnit(annotations,"A0bar","m2 s-1");
        end

        function newFilesUseCanonicalUnitsAndLegacyUnitsStillRestore(testCase)
            path = fullfile(testCase.temporaryFolder,"units.nc");
            wvt = TestPhysicalUnitMetadata.boussinesqTransform();
            ncfile = wvt.writeToFile(char(path),shouldOverwriteExisting=true);
            cleanup = onCleanup(@()TestPhysicalUnitMetadata.closeIfOpen(ncfile));

            testCase.verifyNetCDFUnit(ncfile,"Ap","m s-1");
            testCase.verifyNetCDFUnit(ncfile,"A0","m2 s-1");
            testCase.verifyNetCDFUnit(ncfile,"K2unique","rad2 m-2");
            testCase.verifyNetCDFUnit(ncfile,"iK2unique","1");
            testCase.verifyNetCDFUnit(ncfile,"rho0","kg m-3");

            ncfile.variableWithName("Ap").addAttribute("units","m/s");
            ncfile.variableWithName("A0").addAttribute("units","m^2/s");
            ncfile.variableWithName("K2unique").addAttribute("units","rad^2 m^{-2}");
            ncfile.close();
            clear cleanup

            restored = testCase.verifyWarningFree(@()WVTransform.waveVortexTransformFromFile(char(path)));
            testCase.verifyClass(restored,class(wvt));
            testCase.verifyEqual(real(restored.Ap),real(wvt.Ap));
            testCase.verifyEqual(imag(restored.Ap),imag(wvt.Ap));
            testCase.verifyEqual(real(restored.Am),real(wvt.Am));
            testCase.verifyEqual(imag(restored.Am),imag(wvt.Am));
            testCase.verifyEqual(real(restored.A0),real(wvt.A0));
            testCase.verifyEqual(imag(restored.A0),imag(wvt.A0));
        end

        function generatedDocumentationUsesFormattedCanonicalUnits(testCase)
            boussinesqRoot = fullfile(testCase.repositoryRoot,"docs","classes","transforms","wvtransformboussinesq");
            testCase.verifySubstring(fileread(fullfile(boussinesqRoot,"k2.md")),"units of $$\mathrm{rad^{2}\,m^{-2}}$$.");
            testCase.verifySubstring(fileread(fullfile(boussinesqRoot,"a0t.md")),"units of $$\mathrm{m^{2}\,s^{-1}}$$.");
            testCase.verifySubstring(fileread(fullfile(boussinesqRoot,"p.md")),"units of $$\mathrm{kg\,m^{-1}\,s^{-2}}$$.");
            testCase.verifySubstring(fileread(fullfile(boussinesqRoot,"ap.md")),"units of $$\mathrm{m\,s^{-1}}$$.");
        end
    end

    methods (Access=private)
        function verifyCanonicalAnnotations(testCase,annotations,context)
            allowedUnits = ["1" "degrees_north" "kg m-1 s-2" "kg m-3" "m" "m s-1" "m s-2" "m-1" "m-2" "m-3" "m2" "m2 s-1" "m2 s-2" "m3 s-2" "rad m-1" "rad m-1 s-1" "rad s-1" "rad2 m-2" "rad2 s-2" "s" "s m-1" "s-1" "s-2"];
            allowedEmptyNames = ["A0_Psi_factor" "absTolerance"];
            for iAnnotation = 1:numel(annotations)
                annotation = annotations(iAnnotation);
                if ~isprop(annotation,"units")
                    continue
                end
                units = string(annotation.units);
                name = string(annotation.name);
                if units == ""
                    testCase.verifyTrue(any(name == allowedEmptyNames),context + ": unexpected empty units for " + name);
                else
                    testCase.verifyTrue(any(units == allowedUnits),context + ": unsupported units '" + units + "' for " + name);
                end
            end
        end

        function verifyAnnotationUnit(testCase,annotations,name,expected)
            names = string({annotations.name});
            matches = find(names == name);
            testCase.assertNotEmpty(matches,"Missing annotation " + name);
            testCase.verifyEqual(unique(string({annotations(matches).units})),expected);
        end

        function verifyNetCDFUnit(testCase,ncfile,name,expected)
            variable = ncfile.variableWithName(name);
            testCase.verifyTrue(variable.attributes.isKey("units"),"Missing NetCDF units for " + name);
            testCase.verifyEqual(string(variable.attributes("units")),expected);
        end
    end

    methods (Static,Access=private)
        function transforms = stableTransforms()
            Lxyz = [4000 3000 1000];
            Nxyz = [8 6 5];
            N2 = @(z) 2e-5*exp(z/4000);
            transforms = { ...
                WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false), ...
                WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false), ...
                WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false), ...
                TestPhysicalUnitMetadata.boussinesqTransform(), ...
                WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false), ...
                WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=false)};
        end

        function wvt = boussinesqTransform()
            wvt = WVTransformBoussinesq([4000 3000 1000],[8 6 5],N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
