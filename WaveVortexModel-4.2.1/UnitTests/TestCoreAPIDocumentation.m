classdef TestCoreAPIDocumentation < matlab.unittest.TestCase
    properties
        repositoryRoot (1,1) string
    end

    methods (TestMethodSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test,TestTags="full")
        function abstractTransformListsCurrentFamilies(testCase)
            page = testCase.generatedTransformPage("wvtransform","index.md");
            transformFamilies = [
                "WVTransformConstantStratification"
                "WVTransformHydrostatic"
                "WVTransformBoussinesq"
                "WVTransformStratifiedQG"
                "WVTransformBarotropicQG"
                ];
            for family = transformFamilies'
                testCase.verifySubstring(page,"`" + family + "`");
            end
            testCase.verifyFalse(contains(page,"WVTransformSingleMode"));
            testCase.verifyFalse(contains(page,"Cartesian2DBarotropic"));
        end

        function constructorPagesUseCurrentDeclarations(testCase)
            expectations = {
                "wvtransformconstantstratification", "WVTransformConstantStratification", ["options.N0" "options.isHydrostatic"]
                "wvtransformhydrostatic", "WVTransformHydrostatic", ["options.N2Function" "options.rhoFunction"]
                "wvtransformboussinesq", "WVTransformBoussinesq", ["options.N2Function" "options.rhoFunction"]
                "wvtransformstratifiedqg", "WVTransformStratifiedQG", ["options.N2Function" "options.rhoFunction"]
                "wvtransformbarotropicqg", "WVTransformBarotropicQG", "options.h"
                };
            for iExpectation = 1:size(expectations,1)
                folder = expectations{iExpectation,1};
                className = expectations{iExpectation,2};
                page = testCase.generatedTransformPage(folder,folder + ".md");
                testCase.verifySubstring(page,"wvt = " + className + "(");
                testCase.verifySubstring(page,"new `" + className + "` instance");
                for option = expectations{iExpectation,3}
                    testCase.verifySubstring(page,"`" + option + "`");
                end
                testCase.verifyFalse(contains(page,"WaveVortexTranform"));
                if className ~= "WVTransformHydrostatic"
                    testCase.verifyFalse(contains(page,"WVTransformHydrostatic instance"), ...
                        "A copied Hydrostatic return description remains in " + className + ".");
                end
            end
        end

        function waveVortexCoefficientsAreProminent(testCase)
            abstractPage = testCase.generatedTransformPage("wvtransform","index.md");
            coefficientSection = extractBetween(abstractPage,"+ Inspect wave-vortex coefficients","+ Set and inspect time");
            for name = ["Ap" "Am" "A0"]
                testCase.verifySubstring(coefficientSection,"[`" + name + "`]");
            end

            waveBearing = [
                "wvtransformconstantstratification"
                "wvtransformhydrostatic"
                "wvtransformboussinesq"
                ];
            for folder = waveBearing'
                page = testCase.generatedTransformPage(folder,"index.md");
                for name = ["Ap" "Am" "A0"]
                    testCase.verifySubstring(page,"wvtransform/" + lower(name) + ".html");
                end
                for name = ["Apt" "Amt" "A0t"]
                    testCase.verifySubstring(page,"[`" + name + "`]");
                end
            end

            qgTransforms = ["wvtransformstratifiedqg" "wvtransformbarotropicqg"];
            for folder = qgTransforms
                page = testCase.generatedTransformPage(folder,"index.md");
                testCase.verifySubstring(page,"wvtransform/a0.html");
                testCase.verifySubstring(page,"no active `Ap`, `Am`, `Apt`, or `Amt` content");
                testCase.verifySubstring(page,"[`A0t`]");
            end
        end

        function lowLevelCoefficientArraysAreDeveloperReference(testCase)
            page = testCase.generatedTransformPage("wvtransformconstantstratification","index.md");
            developerSection = extractAfter(page,"## Developer Topics");
            testCase.verifySubstring(developerSection,"Projection and reconstruction coefficients");
            for name = ["A0N" "A0U" "A0V" "UAp" "VAm" "NA0"]
                testCase.verifySubstring(developerSection,"[`" + name + "`]");
            end

            source = fileread(fullfile(testCase.repositoryRoot,"@WVTransform","detailedDescriptions","ApU.md"));
            testCase.verifySubstring(source,"projection coefficients");
            testCase.verifyFalse(contains(source,"sorting matrix"));
        end

        function reviewedPagesContainNoObsoleteReferenceText(testCase)
            relativePaths = [
                "classes/transforms/wvtransform/index.md"
                "classes/transforms/wvtransformconstantstratification/index.md"
                "classes/transforms/wvtransformhydrostatic/index.md"
                "classes/transforms/wvtransformboussinesq/index.md"
                "classes/transforms/wvtransformstratifiedqg/index.md"
                "classes/transforms/wvtransformbarotropicqg/index.md"
                "classes/wvmodel/index.md"
                ];
            forbidden = [
                "WVTransformSingleMode"
                "Cartesian2DBarotropic"
                "WaveVortexTranform"
                "sorting matrix"
                ];
            for relativePath = relativePaths'
                page = fileread(fullfile(testCase.repositoryRoot,"docs",relativePath));
                for phrase = forbidden'
                    testCase.verifyFalse(contains(page,phrase),relativePath + " contains " + phrase + ".");
                end
            end

            interpolationPage = testCase.generatedTransformPage("wvtransform","variableatpositionwithname.md");
            testCase.verifySubstring(interpolationPage,"`linear`");
            testCase.verifySubstring(interpolationPage,"`spline`");
            testCase.verifyFalse(contains(interpolationPage,"`exact`"));
            testCase.verifyFalse(contains(interpolationPage,"`finufft`"));
        end

        function generatedReferenceUsesIntentionalTopics(testCase)
            indexPages = dir(fullfile(testCase.repositoryRoot,"docs","classes","**","index.md"));
            testCase.verifyNotEmpty(indexPages);
            vagueTopics = [
                newline + "+ Other" + newline
                newline + "+ Internal" + newline
                newline + "+ Utility function" + newline
                newline + "+ Properties" + newline
                newline + "+ Index gymnastics" + newline
                ];
            for iPage = 1:numel(indexPages)
                pagePath = fullfile(indexPages(iPage).folder,indexPages(iPage).name);
                page = string(fileread(pagePath));
                for topic = vagueTopics'
                    testCase.verifyFalse(contains(page,topic),pagePath + " contains the vague topic " + strtrim(topic) + ".");
                end
            end
        end

        function landingPagesLeadWithUserTasks(testCase)
            transformPage = testCase.generatedTransformPage("wvtransform","index.md");
            modelPage = string(fileread(fullfile(testCase.repositoryRoot,"docs","classes","wvmodel","index.md")));
            testCase.verifyFalse(contains(transformPage,"Stable"));
            testCase.verifyFalse(any(contains(modelPage,["Stable","Experimental","adaptive-cell"])));
            testCase.verifySubstring(modelPage,"`setupIntegrator` to change time-stepping settings");

            versionPage = testCase.generatedTransformPage("wvtransform","version.md");
            testCase.verifySubstring(versionPage,"Installed WaveVortexModel version");
            testCase.verifyFalse(contains(versionPage,"mpackage.json"));
            packageSection = extractAfter(transformPage,"+ Get package information");
            testCase.verifySubstring(packageSection,"[`version`]");
        end

        function capabilityPageAvoidsReleaseStatusJargon(testCase)
            canonicalPath = fullfile(testCase.repositoryRoot,"Documentation","WebsiteDocumentation", ...
                "users-guide","supported-features.md");
            page = string(fileread(canonicalPath));
            testCase.verifySubstring(page,"# Capabilities and limitations");
            testCase.verifyFalse(contains(page,"Status definitions"));
            testCase.verifyFalse(contains(page,"| Status |"));
            testCase.verifyFalse(contains(page,"Stable"));
        end

        function documentedConstructorsExecute(testCase)
            Lxyz = [40e3 30e3 2e3];
            Nxyz = [8 6 5];
            N2 = @(z)(5.2e-3)^2*ones(size(z));

            transforms = {
                WVTransformConstantStratification(Lxyz,Nxyz,N0=5.2e-3,latitude=45)
                WVTransformConstantStratification(Lxyz,Nxyz,N0=5.2e-3,latitude=45,isHydrostatic=true)
                WVTransformHydrostatic(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformBoussinesq(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformStratifiedQG(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),h=0.8,latitude=45)
                };
            expectedClasses = [
                "WVTransformConstantStratification"
                "WVTransformConstantStratification"
                "WVTransformHydrostatic"
                "WVTransformBoussinesq"
                "WVTransformStratifiedQG"
                "WVTransformBarotropicQG"
                ];
            testCase.verifyEqual(string(cellfun(@class,transforms,UniformOutput=false)),expectedClasses);

            linearModel = WVModel(transforms{1},shouldUseLinearDynamics=true);
            transforms{2}.addForcing(WVAdaptiveDamping(transforms{2}));
            nonlinearModel = WVModel(transforms{2});
            testCase.verifyTrue(linearModel.isDynamicsLinear);
            testCase.verifyFalse(nonlinearModel.isDynamicsLinear);
        end

        function coefficientClaimsMatchTransformState(testCase)
            wvt = WVTransformConstantStratification([40e3 30e3 2e3],[8 6 5],N0=5.2e-3,latitude=45);
            testCase.verifyEqual(size(wvt.Ap),size(wvt.Am));
            testCase.verifyEqual(size(wvt.Ap),size(wvt.A0));
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("Ap").units),"m/s");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("Am").units),"m/s");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("A0").units),"m^2 s^{-1}");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("A0t").units),"m^2 s^{-1}");

            waveIndex = find(wvt.waveComponent.maskAp,1);
            wvt.Ap(waveIndex) = 1.25 - 0.5i;
            wvt.t0 = 7;
            wvt.t = 19;
            expected = wvt.Ap(waveIndex)*exp(1i*wvt.Omega(waveIndex)*(wvt.t-wvt.t0));
            testCase.verifyEqual(wvt.Apt(waveIndex),expected,AbsTol=10*eps(abs(expected)));
            testCase.verifyEqual(wvt.Ap(waveIndex),1.25 - 0.5i);

            wvt.initWithInertialMotions(@(z)0.1*ones(size(z)),@(z)zeros(size(z)));
            maskAp = logical(wvt.inertialComponent.maskAp);
            maskAm = logical(wvt.inertialComponent.maskAm);
            testCase.verifyEqual(wvt.Am(maskAm),conj(wvt.Ap(maskAp)));

            qg = WVTransformBarotropicQG([40e3 30e3],[8 6],h=0.8,latitude=45);
            testCase.verifyFalse(qg.hasVariableWithName("Ap"));
            testCase.verifyFalse(qg.hasVariableWithName("Am"));
            testCase.verifyTrue(qg.hasVariableWithName("A0"));
            testCase.verifyEqual(qg.A0t,qg.A0);
        end

        function authoredScientificContextIsPreserved(testCase)
            coefficientFiles = ["Ap.md" "Am.md" "A0.md"];
            for filename = coefficientFiles
                source = testCase.coefficientSidecar(filename);
                testCase.verifySubstring(source,"| Vertical mode | $$K_h=0$$ | $$K_h>0$$ |");
                testCase.verifySubstring(source,"10.1017/jfm.2020.995");
                testCase.verifySubstring(source,"10.48550/arXiv.2403.20269");
                testCase.verifySubstring(source,"executable definitions");

                generatedPage = testCase.generatedTransformPage("wvtransform",lower(erase(filename,".md")) + ".md");
                testCase.verifySubstring(generatedPage,"| Vertical mode | $$K_h=0$$ | $$K_h>0$$ |");
                testCase.verifySubstring(generatedPage,"10.1017/jfm.2020.995");
            end
            testCase.verifySubstring(testCase.coefficientSidecar("Ap.md"),"internal gravity wave, $$+\omega$$");
            testCase.verifySubstring(testCase.coefficientSidecar("Am.md"),"internal gravity wave, $$-\omega$$");
            a0 = testCase.coefficientSidecar("A0.md");
            testCase.verifySubstring(a0,"barotropic geostrophic");
            testCase.verifySubstring(a0,"baroclinic geostrophic");
            testCase.verifySubstring(a0,"mean density anomaly");
            testCase.verifySubstring(a0,"$$K_h=0,j=0$$ location carries neither");
        end

        function authoredDesignRationaleIsPreserved(testCase)
            expectations = {
                fullfile("classes","model-output","wvmodeloutputgroup","index.md"), "resolve the buoyancy frequency"
                fullfile("classes","model-output","wvmodeloutputfile","index.md"), "writes nothing until a group is"
                fullfile("classes","forcing","wvforcing","index.md"), "Forcing is applied in three stages"
                fullfile("classes","flow-components","wvflowcomponent","index.md"), "analytical mode, identified"
                fullfile("classes","operations-and-annotations","wvvariableannotation","index.md"), "keeps equations and tables out"
                fullfile("classes","transforms","wvtransform","index.md"), "No temporal filter is"
                fullfile("classes","wvmodel","index.md"), "Model output is assembled in three layers"
                };
            for iExpectation = 1:size(expectations,1)
                page = string(fileread(fullfile(testCase.repositoryRoot,"docs",expectations{iExpectation,1})));
                testCase.verifySubstring(page,expectations{iExpectation,2});
            end
        end

        function coefficientOccupancyMatchesFlowComponentMasks(testCase)
            transforms = testCase.supportedTransformFixtures();
            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                components = wvt.primaryFlowComponents;
                assignedAp = zeros(wvt.spectralMatrixSize);
                assignedAm = zeros(wvt.spectralMatrixSize);
                assignedA0 = zeros(wvt.spectralMatrixSize);

                for iComponent = 1:numel(components)
                    component = components(iComponent);
                    maskAp = logical(component.maskAp);
                    maskAm = logical(component.maskAm);
                    maskA0 = logical(component.maskA0);
                    assignedAp = assignedAp + maskAp;
                    assignedAm = assignedAm + maskAm;
                    assignedA0 = assignedA0 + maskA0;

                    switch string(component.shortName)
                        case "wave"
                            expected = wvt.Kh > 0 & wvt.J > 0;
                            testCase.verifyEqual(maskAp,expected,class(wvt) + " wave Ap occupancy");
                            testCase.verifyEqual(maskAm,expected,class(wvt) + " wave Am occupancy");
                            testCase.verifyFalse(any(maskA0,"all"));
                        case "inertial"
                            expected = wvt.Kh == 0;
                            testCase.verifyEqual(maskAp,expected,class(wvt) + " inertial Ap occupancy");
                            testCase.verifyEqual(maskAm,expected,class(wvt) + " inertial Am occupancy");
                            testCase.verifyFalse(any(maskA0,"all"));
                        case "geostrophic"
                            expected = wvt.Kh > 0;
                            testCase.verifyEqual(maskA0,expected,class(wvt) + " geostrophic A0 occupancy");
                            testCase.verifyFalse(any(maskAp,"all"));
                            testCase.verifyFalse(any(maskAm,"all"));
                        case "mda"
                            expected = wvt.Kh == 0 & wvt.J > 0;
                            testCase.verifyEqual(maskA0,expected,class(wvt) + " MDA A0 occupancy");
                            testCase.verifyFalse(any(maskAp,"all"));
                            testCase.verifyFalse(any(maskAm,"all"));
                        otherwise
                            testCase.assertFail("Unexpected primary flow component " + string(component.shortName));
                    end
                end

                testCase.verifyLessThanOrEqual(assignedAp,ones(size(assignedAp)),class(wvt) + " overlapping Ap masks");
                testCase.verifyLessThanOrEqual(assignedAm,ones(size(assignedAm)),class(wvt) + " overlapping Am masks");
                testCase.verifyLessThanOrEqual(assignedA0,ones(size(assignedA0)),class(wvt) + " overlapping A0 masks");
                testCase.verifyFalse(any(assignedA0(wvt.Kh == 0 & wvt.J == 0)),class(wvt) + " assigned the empty A0 origin");
            end
        end
    end

    methods (Access=private)
        function page = generatedTransformPage(testCase,folder,file)
            page = string(fileread(fullfile(testCase.repositoryRoot,"docs","classes","transforms",folder,file)));
        end

        function source = coefficientSidecar(testCase,filename)
            source = string(fileread(fullfile(testCase.repositoryRoot,"@WVTransform","detailedDescriptions",filename)));
        end

        function transforms = supportedTransformFixtures(~)
            Lxyz = [40e3 30e3 2e3];
            Nxyz = [8 6 5];
            N2 = @(z)(5.2e-3)^2*ones(size(z));
            transforms = {
                WVTransformConstantStratification(Lxyz,Nxyz,N0=5.2e-3,latitude=45)
                WVTransformConstantStratification(Lxyz,Nxyz,N0=5.2e-3,latitude=45,isHydrostatic=true)
                WVTransformHydrostatic(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformBoussinesq(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformStratifiedQG(Lxyz,Nxyz,N2Function=N2,latitude=45)
                WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),h=0.8,latitude=45)
                };
        end
    end
end
