classdef TestForcingDocumentation < matlab.unittest.TestCase
    properties
        repositoryRoot (1,1) string
    end

    methods (TestMethodSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test, TestTags="full")
        function suppliedForcingPagesDescribeCurrentContracts(testCase)
            expectations = {
                "wvnonlinearadvection", false, "self = WVNonlinearAdvection(wvt)", ["installs this forcing by default" "using linear evolution"]
                "wvfixedamplitudeforcing", false, "self = WVFixedAmplitudeForcing(wvt,options)", ["`Apbar`" "`A0_indices`" "10^{-6}"]
                "wvnarrowbandgeostrophicforcing", false, "self = WVNarrowBandGeostrophicForcing(wvt,options)", ["`modelSpectrum`" "`initialPV`" "consume random numbers"]
                "wvbottomfrictionlinear", false, "self = WVBottomFrictionLinear(wvt,options)", ["`r`" "inverse seconds" "200-day"]
                "wvbottomfrictionquadratic", false, "self = WVBottomFrictionQuadratic(wvt,options)", ["`Cd`" "dimensionless" "4000"]
                "wvadaptivedamping", true, "self = WVAdaptiveDamping(wvt)", ["wvt.uvMax*damp" "k_no_damp" "significant damping"]
                "wvhorizontaldamping", true, "self = WVHorizontalDamping(wvt,options)", ["`nu`" "`kappa`" "square meters per second"]
                "wvverticaldamping", true, "self = WVVerticalDamping(wvt,options)", ["`nu`" "`kappa`" "square meters per second"]
                "wvverticaldiffusivity", true, "self = WVVerticalDiffusivity(wvt,options)", ["`kappa_z`" "`shouldForceMeanDensityAnomaly`" "barotropic QG"]
                "wvantialiasing", true, "self = WVAntialiasing(wvt,options)", ["`Nj`" "shouldAntialias=false" "discarded"]
                "wvthermaldamping", true, "self = WVThermalDamping(wvt,options)", ["`alpha`" "alpha/wvt.Lr2" "Scott and Dritschel"]
                };

            for iExpectation = 1:size(expectations,1)
                folder = expectations{iExpectation,1};
                isClosure = expectations{iExpectation,2};
                constructorDeclaration = expectations{iExpectation,3};
                requiredText = expectations{iExpectation,4};
                overview = testCase.generatedPage(folder,"index.md",isClosure);
                constructor = testCase.generatedPage(folder,folder + ".md",isClosure);
                reference = testCase.generatedClassReference(folder,isClosure);
                testCase.verifySubstring(overview,"### Example",folder + " lacks an example.");
                testCase.verifySubstring(constructor,constructorDeclaration,folder + " has the wrong constructor declaration.");
                for token = requiredText
                    testCase.verifySubstring(reference,token,folder + " lacks " + token + ".");
                end
            end
        end

        function documentationContainsNoKnownCopiedDefects(testCase)
            generatedRoot = fullfile(testCase.repositoryRoot,"docs","classes","forcing");
            pages = dir(fullfile(generatedRoot,"**","*.md"));
            forbidden = [
                "WVNonlinearFlux"
                "WVAdaptiveViscosity"
                "WVVerticalScalarDiffusivity"
                "WaveVortexTranform"
                "comptued"
                "diffusivty"
                "C_dd"
                "WVNonlinearAdvection(wvt,options)"
                ];
            for iPage = 1:numel(pages)
                page = string(fileread(fullfile(pages(iPage).folder,pages(iPage).name)));
                for token = forbidden'
                    testCase.verifyFalse(contains(page,token),pages(iPage).name + " contains obsolete text " + token + ".");
                end
            end
        end

        function publishedExamplesExecute(testCase)
            classes = [
                "WVNonlinearAdvection"
                "WVFixedAmplitudeForcing"
                "WVNarrowBandGeostrophicForcing"
                "WVBottomFrictionLinear"
                "WVBottomFrictionQuadratic"
                "WVAdaptiveDamping"
                "WVHorizontalDamping"
                "WVVerticalDamping"
                "WVVerticalDiffusivity"
                "WVAntialiasing"
                "WVThermalDamping"
                ];
            for className = classes'
                sourcePath = fullfile(testCase.repositoryRoot,"Forcing",className + ".m");
                code = testCase.exampleCodeFromSource(sourcePath);
                testCase.assertNotEmpty(code,className + " has no executable example.");
                testCase.verifyWarningFree(@()testCase.evaluateExample(code),className + " example emitted a warning.");
            end
        end

        function documentedDerivedQuantitiesMatchImplementation(testCase)
            Lxyz = [40e3 30e3 2e3];
            Nxyz = [8 6 5];
            wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=5.2e-3,latitude=45,isHydrostatic=true,shouldAntialias=false);

            linear = WVBottomFrictionLinear(wvt,r=2e-7);
            testCase.verifyEqual(linear.r_scaled,linear.r*wvt.Lz/wvt.z_int(1),RelTol=10*eps);

            quadratic = WVBottomFrictionQuadratic(wvt,Cd=2e-3);
            testCase.verifyEqual(quadratic.cd,quadratic.Cd/wvt.z_int(1),RelTol=10*eps);

            horizontal = WVHorizontalDamping(wvt);
            testCase.verifyEqual(horizontal.nu,1e-4);
            testCase.verifyEqual(horizontal.kappa,1e-6);

            vertical = WVVerticalDamping(wvt);
            testCase.verifyEqual(vertical.nu,5e-4);
            testCase.verifyEqual(vertical.kappa,1e-6);

            diffusivity = WVVerticalDiffusivity(wvt);
            testCase.verifyEqual(diffusivity.kappa_z,1e-5);
            testCase.verifyTrue(diffusivity.shouldForceMeanDensityAnomaly);

            adaptive = WVAdaptiveDamping(wvt);
            testCase.verifySize(adaptive.damp,wvt.spectralMatrixSize);
            testCase.verifyLessThan(adaptive.k_no_damp,adaptive.k_damp);
            testCase.verifyLessThan(adaptive.j_no_damp,adaptive.j_damp);
            testCase.verifyEqual(adaptive.dampingTimeScale(),1/max(abs(adaptive.damp),[],"all"),RelTol=10*eps);

            antialias = WVAntialiasing(wvt);
            testCase.verifyEqual(antialias.Nj,floor(2*wvt.Nj/3));
            testCase.verifySize(antialias.M,wvt.spectralMatrixSize);
            testCase.verifyEqual(antialias.effectiveJMax(),antialias.Nj-1);

            qg = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),h=0.8,latitude=45,shouldAntialias=false);
            qgLinear = WVBottomFrictionLinear(qg,r=2e-7);
            qgQuadratic = WVBottomFrictionQuadratic(qg,Cd=2e-3);
            thermal = WVThermalDamping(qg,alpha=3e-8);
            testCase.verifyEqual(qgLinear.r_scaled,qgLinear.r);
            testCase.verifyEqual(qgQuadratic.cd,qgQuadratic.Cd/4000,RelTol=10*eps);
            testCase.verifyEqual(thermal.alpha_scaled,thermal.alpha./qg.Lr2,RelTol=10*eps);
        end

        function fixedAmplitudeHelperDocumentsCurrentBehavior(testCase)
            page = testCase.generatedPage("wvfixedamplitudeforcing","setnarrowbandgeostrophicforcing.md",false);
            required = [
                "`k_f`"
                "`j_f`"
                "`k_r`"
                "`u_rms`"
                "`initialPV`"
                "`r`"
                "narrow-band"
                "full-spectrum"
                "WVNarrowBandGeostrophicForcing"
                "Deprecated 4.x helper"
                ];
            for token = required'
                testCase.verifySubstring(page,token);
            end

            wvt = WVTransformConstantStratification([40e3 30e3 2e3],[8 6 5],N0=5.2e-3,latitude=45,isHydrostatic=true);
            forcing = WVFixedAmplitudeForcing(wvt,name="selected-geostrophic-mode");
            A0bar = zeros(wvt.spectralMatrixSize);
            selectedIndex = find(wvt.geostrophicComponent.maskA0,1);
            A0bar(selectedIndex) = 1;
            forcing.setGeostrophicForcingCoefficients(A0bar);
            testCase.verifyEqual(forcing.A0_indices,uint64(selectedIndex));
            testCase.verifyEqual(forcing.A0bar,1);
        end
    end

    methods (Access=private)
        function page = generatedPage(testCase,folder,filename,isClosure)
            root = fullfile(testCase.repositoryRoot,"docs","classes","forcing");
            if isClosure
                root = fullfile(root,"closures");
            end
            page = string(fileread(fullfile(root,folder,filename)));
        end

        function reference = generatedClassReference(testCase,folder,isClosure)
            root = fullfile(testCase.repositoryRoot,"docs","classes","forcing");
            if isClosure
                root = fullfile(root,"closures");
            end
            pages = dir(fullfile(root,folder,"*.md"));
            reference = "";
            for iPage = 1:numel(pages)
                reference = reference + newline + string(fileread(fullfile(pages(iPage).folder,pages(iPage).name)));
            end
        end

        function code = exampleCodeFromSource(~,sourcePath)
            source = fileread(sourcePath);
            match = regexp(source,'(?s)% ### Example.*?% ```matlab\s*(?<code>.*?)\s*% ```','names','once');
            if isempty(match)
                code = "";
                return
            end
            code = string(regexprep(match.code,'(?m)^\s*% ?',''));
        end

        function evaluateExample(~,code)
            evalc(char(code));
        end
    end
end
