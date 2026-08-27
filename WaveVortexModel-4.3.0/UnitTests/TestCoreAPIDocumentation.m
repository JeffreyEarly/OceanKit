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
            coefficientSection = extractBetween(abstractPage,"+ Inspect wave-vortex coefficients","+ Create a related transform");
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

        function concreteTransformsUseTaskOrientedTopics(testCase)
            expected = [
                "Create and restore a transform"
                "Inspect the domain"
                "Initialize the flow"
                "Evaluate physical fields"
                "Manage forcing and closures"
                "Analyze the flow"
                "Analyze energy"
                "Save transform state"
                "Convert representations"
                "Differentiate and integrate fields"
                "Inspect flow components"
                "Inspect wave-vortex coefficients"
                "Create a related transform"
                "Extend a transform"
                "Get package information"
                ];
            folders = [
                "wvtransformconstantstratification"
                "wvtransformhydrostatic"
                "wvtransformboussinesq"
                "wvtransformstratifiedqg"
                "wvtransformbarotropicqg"
                ];
            for folder = folders'
                page = testCase.generatedTransformPage(folder,"index.md");
                testCase.verifyEqual(testCase.topLevelTopics(page),expected,folder);
                testCase.verifyLessThan(strfind(page,"+ Initialize the flow"), ...
                    strfind(page,"+ Inspect wave-vortex coefficients"),folder);
                testCase.verifyNoEmptyTopicBranches(page,folder);
            end
        end

        function transformDomainTopicsHaveScientificOrder(testCase)
            page = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            domain = extractBetween(page,"+ Inspect the domain","+ Initialize the flow");
            testCase.verifyTextOrder(domain,["[`x`]" "[`y`]" "[`z`]"]);
            testCase.verifyTextOrder(domain,["[`X`]" "[`Y`]" "[`Z`]" "[`xyzGrid`]"]);
            testCase.verifyTextOrder(domain,["[`z_int`]" "[`volumeIntegral`]"]);
            testCase.verifyTextOrder(domain,["[`k`]" "[`l`]" "[`j`]"]);
            testCase.verifyTextOrder(domain,["[`K`]" "[`L`]" "[`J`]" "[`kljGrid`]"]);
            testCase.verifyTextOrder(domain,["[`dk`]" "[`dl`]"]);
            testCase.verifyTextOrder(domain,["[`Kh`]" "[`K2`]"]);
            testCase.verifyTextOrder(domain,["[`FMatrix`]" "[`FinvMatrix`]" "[`GMatrix`]" "[`GinvMatrix`]"]);
            testCase.verifyFalse(any(contains(domain,["[`kAxis`]" "[`lAxis`]" "[`transformToKLAxes`]" ])));
            for name = ["f" "g" "dLnN2" "rho0" "planetaryRadius" "rotationRate"]
                testCase.verifySubstring(domain,"[`" + name + "`]");
            end

            barotropic = testCase.generatedTransformPage("wvtransformbarotropicqg","index.md");
            barotropicDomain = extractBetween(barotropic,"+ Inspect the domain","+ Initialize the flow");
            testCase.verifyTextOrder(barotropicDomain,["[`k`]" "[`l`]"]);
            testCase.verifyTextOrder(barotropicDomain,["[`K`]" "[`L`]" "[`klGrid`]"]);

            analysis = extractBetween(page,"+ Analyze the flow","+ Analyze energy");
            testCase.verifyTextOrder(analysis,["[`kAxis`]" "[`lAxis`]" "[`transformToKLAxes`]"]);
            testCase.verifyTextOrder(analysis,["[`kRadial`]" "[`transformToRadialWavenumber`]"]);
            testCase.verifyFalse(contains(analysis,"Pseudo-radial wavenumber"));
            testCase.verifyFalse(contains(analysis,"transformToOmegaAxis"));
        end

        function transformTopicsGroupForcingEnergyAndComponents(testCase)
            page = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            forcing = extractBetween(page,"+ Manage forcing and closures","+ Analyze the flow");
            testCase.verifyTextOrder(forcing,["Configure forcing" "Inspect forcing and closures" "Summarize forcing"]);
            testCase.verifySubstring(forcing,"[`addForcing`]");
            testCase.verifySubstring(forcing,"[`forcingNames`]");
            testCase.verifySubstring(forcing,"[`summarizeForcing`]");

            energy = extractBetween(page,"+ Analyze energy","+ Save transform state");
            testCase.verifyTextOrder(energy,["Component energy" "Total energy" "Energy summaries"]);
            testCase.verifySubstring(energy,"[`totalEnergyOfFlowComponent`]");
            testCase.verifySubstring(energy,"[`totalEnergy`]");
            testCase.verifySubstring(energy,"[`summarizeEnergyContent`]");

            components = extractBetween(page,"+ Inspect flow components","+ Inspect wave-vortex coefficients");
            testCase.verifyTextOrder(components,["Primary flow components" "Registered and combined components" "Summarize flow components"]);
            testCase.verifyTextOrder(components,["[`waveComponent`]" "[`inertialComponent`]" "[`geostrophicComponent`]" "[`mdaComponent`]"]);

            for folder = ["wvtransformstratifiedqg" "wvtransformbarotropicqg"]
                qgPage = testCase.generatedTransformPage(folder,"index.md");
                qgComponents = extractBetween(qgPage,"+ Inspect flow components","+ Inspect wave-vortex coefficients");
                testCase.verifySubstring(qgComponents,"[`geostrophicComponent`]");
                testCase.verifyFalse(any(contains(qgComponents,["[`waveComponent`]" "[`inertialComponent`]" "[`mdaComponent`]" ])));
            end
        end

        function transformDeveloperGeometryHasIntentionalSubtopics(testCase)
            page = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            developer = extractAfter(page,"## Developer Topics");
            geometry = extractBetween(developer,"+ Geometry and mode indexing","+ Spectral transforms and operators");
            expectedSubtopics = [
                "Mode numbers and validity"
                "Linear-index conversion"
                "DFT and WV layout metadata"
                "Layout conversion"
                "Masks and Hermitian bookkeeping"
                ];
            testCase.verifyTextOrder(geometry,expectedSubtopics);

            movedMembers = ["Nk_dft" "Nl_dft" "k_dft" "l_dft" "kl" ...
                "dftConjugateIndices2D" "dftPrimaryIndices2D" ...
                "indicesOfFourierConjugates" "shouldExcludeConjugates" ...
                "shouldExcludeNyquist" "isHermitian" "setConjugateToUnity"];
            for name = movedMembers
                testCase.verifySubstring(geometry,"[`" + name + "`]");
            end

            classInternals = extractAfter(developer,"+ Class internals");
            for name = movedMembers
                testCase.verifyFalse(contains(classInternals,"[`" + name + "`]"));
            end
        end

        function transformSummariesAreConciseAndSpecific(testCase)
            page = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            expectedSummaries = {
                "Ap", "Positive-frequency wave–vortex coefficient array."
                "Am", "Negative-frequency wave–vortex coefficient array."
                "A0", "Zero-frequency wave–vortex coefficient array."
                "k", "Compact `Nkl`-by-1 x-wavenumber vector in $$\mathrm{rad\,m^{-1}}$$."
                "l", "Compact `Nkl`-by-1 y-wavenumber vector in $$\mathrm{rad\,m^{-1}}$$."
                "j", "Dimensionless `Nj`-by-1 vertical-mode index vector."
                "kAxis", "Centered `Nx`-by-1 x-wavenumber axis in $$\mathrm{rad\,m^{-1}}$$."
                "lAxis", "Centered `Ny`-by-1 y-wavenumber axis in $$\mathrm{rad\,m^{-1}}$$."
                "rho_bar", "Current horizontally averaged density, `[Nz 1]`, in $$\mathrm{kg\,m^{-3}}$$."
                "rho_nm", "Diagnosed no-motion density profile, `[Nz 1]`, in $$\mathrm{kg\,m^{-3}}$$."
                "rho_nm0", "Reference no-motion density profile, `[Nz 1]`, in $$\mathrm{kg\,m^{-3}}$$."
                "FMatrix", "Transformation matrix $$F$$ projecting F-grid values onto vertical modes; shape `[Nj Nz]`."
                "FinvMatrix", "Transformation matrix $$F^{-1}$$ reconstructing F-grid values from vertical modes; shape `[Nz Nj]`."
                "GMatrix", "Transformation matrix $$G$$ projecting G-grid values onto vertical modes; shape `[Nj Nz]`."
                "GinvMatrix", "Transformation matrix $$G^{-1}$$ reconstructing G-grid values from vertical modes; shape `[Nz Nj]`."
                };
            for iExpectation = 1:size(expectedSummaries,1)
                testCase.verifyEqual(testCase.memberSummary(page,expectedSummaries{iExpectation,1}),expectedSummaries{iExpectation,2});
            end

            coefficientNames = ["A0N" "A0U" "A0V" "A0Z" "ApmD" "ApmN" ...
                "NA0" "NAm" "NAp" "PA0" "UA0" "UAm" "UAp" ...
                "VA0" "VAm" "VAp" "WAm" "WAp"];
            for name = coefficientNames
                summary = testCase.memberSummary(page,name);
                testCase.verifyTrue(startsWith(summary,["Projects" "Reconstructs"]),name + " has an inconsistent summary: " + summary);
                testCase.verifyLessThanOrEqual(strlength(summary),80,name + " has an overly long topic summary.");
            end
        end

        function staleTransformLevelSpectralBinningIsAbsent(testCase)
            folders = [
                "wvtransformconstantstratification"
                "wvtransformhydrostatic"
                "wvtransformboussinesq"
                "wvtransformstratifiedqg"
                ];
            retiredMembers = [
                "kPseudoRadial"
                "transformToPseudoRadialWavenumber"
                "transformToPseudoRadialWavenumberA0"
                "transformToPseudoRadialWavenumberApm"
                "transformToOmegaAxis"
                ];
            for folder = folders'
                page = testCase.generatedTransformPage(folder,"index.md");
                for member = retiredMembers'
                    testCase.verifyFalse(contains(page,"[`" + member + "`]"),folder + " still documents " + member + ".");
                end
                for member = ["kRadial" "transformToRadialWavenumber"]
                    testCase.verifySubstring(page,"[`" + member + "`]",folder + " lost " + member + ".");
                end
            end

            retiredPageNames = [
                "kpseudoradial.md"
                "transformtopseudoradialwavenumber.md"
                "transformtopseudoradialwavenumbera0.md"
                "transformtopseudoradialwavenumberapm.md"
                "transformtoomegaaxis.md"
                ];
            generatedRoot = fullfile(testCase.repositoryRoot,"docs","classes","transforms");
            for folder = ["wvtransformhydrostatic" "wvtransformboussinesq" "wvtransformstratifiedqg"]
                for pageName = retiredPageNames'
                    testCase.verifyFalse(isfile(fullfile(generatedRoot,folder,pageName)),folder + " still contains " + pageName + ".");
                end
            end

            N2 = @(z)(5.2e-3)^2*ones(size(z));
            transforms = {
                WVTransformConstantStratification([4000 3000 1000],[8 6 5],N0=5.2e-3,latitude=45)
                WVTransformHydrostatic([4000 3000 1000],[8 6 5],N2Function=N2,latitude=45)
                WVTransformBoussinesq([4000 3000 1000],[8 6 5],N2Function=N2,latitude=45)
                WVTransformStratifiedQG([4000 3000 1000],[8 6 5],N2Function=N2,latitude=45)
                };
            for transform = transforms'
                wvt = transform{1};
                for member = retiredMembers'
                    testCase.verifyFalse(isprop(wvt,member),class(wvt) + " still has property " + member + ".");
                    testCase.verifyFalse(ismethod(wvt,member),class(wvt) + " still has method " + member + ".");
                end
                radial = wvt.transformToRadialWavenumber(ones(wvt.spectralMatrixSize));
                testCase.verifyEqual(size(radial),[wvt.Nj length(wvt.kRadial)]);
                testCase.verifyTrue(all(isfinite(radial),"all"));
            end

            for transform = transforms(1:3)'
                wvt = transform{1};
                phase = wvt.variableWithName("phase");
                conjPhase = wvt.variableWithName("conjPhase");
                testCase.verifySize(wvt.Omega,wvt.spectralMatrixSize);
                testCase.verifySize(wvt.iOmega,wvt.spectralMatrixSize);
                testCase.verifyTrue(all(isfinite(wvt.Omega),"all"));
                testCase.verifyTrue(all(isfinite(wvt.iOmega),"all"));
                testCase.verifyEqual(wvt.iOmega,1i*wvt.Omega);
                testCase.verifyEqual(conjPhase,conj(phase));
            end
        end

        function concreteTransformsShowUsefulInheritedSurface(testCase)
            hydrostatic = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            for name = ["variableWithName" "addForcing" "writeToFile" "Ap" "flowComponents"]
                testCase.verifySubstring(hydrostatic, ...
                    "/classes/transforms/wvtransformhydrostatic/" + lower(name) + ".html");
            end
            testCase.verifyFalse(contains(hydrostatic,"geometryfromfile.html"));
            testCase.verifySubstring(extractAfter(hydrostatic,"## Developer Topics"),"shouldExcludeConjugates");
            testCase.verifySubstring(extractAfter(hydrostatic,"## Developer Topics"),"enstrophyFluxFromF0");

            for folder = ["wvtransformstratifiedqg" "wvtransformbarotropicqg"]
                page = testCase.generatedTransformPage(folder,"index.md");
                testCase.verifyFalse(contains(page,"/" + folder + "/ap.html"));
                testCase.verifyFalse(contains(page,"/" + folder + "/am.html"));
                testCase.verifySubstring(page,"/" + folder + "/a0.html) Zero-frequency geostrophic coefficients.");
            end

            barotropic = testCase.generatedTransformPage("wvtransformbarotropicqg","index.md");
            for name = ["z" "z_" "lz" "nz" "j" "j_" "nj" "kljgrid" "effectivejmax" ...
                    "initwithuveta" "transformuvetatowavevortex"]
                testCase.verifyFalse(contains(barotropic, ...
                    "/classes/transforms/wvtransformbarotropicqg/" + name + ".html"));
            end
            testCase.verifySubstring(barotropic,"Equivalent depth and deformation scale");
            testCase.verifySubstring(barotropic,"[`psi`]");
        end

        function concreteTransformMembersHaveMeaningfulSummaries(testCase)
            folders = [
                "wvtransformconstantstratification"
                "wvtransformhydrostatic"
                "wvtransformboussinesq"
                "wvtransformstratifiedqg"
                "wvtransformbarotropicqg"
                ];
            placeholders = ["dimension","[Nj 1]","This preserves"];
            for folder = folders'
                page = testCase.generatedTransformPage(folder,"index.md");
                topicBlock = extractBetween(page,"## Topics","## Developer Topics");
                memberLines = splitlines(topicBlock);
                memberLines = memberLines(contains(memberLines,"+ [`"));
                testCase.assertNotEmpty(memberLines,folder + " has no documented members.");
                for iMember = 1:numel(memberLines)
                    name = string(regexp(memberLines(iMember),'\[`(?<name>[^`]+)`\]','names','once').name);
                    summary = strtrim(regexprep(memberLines(iMember),'^.*\)\s*',''));
                    testCase.verifyNotEmpty(summary,folder + "." + name + " has no summary.");
                    testCase.verifyFalse(any(summary == placeholders),folder + "." + name + " has placeholder documentation: " + summary);
                    testCase.verifyFalse(endsWith(summary,[" the"," a"," an"," and"," or"," with"," when"]),folder + "." + name + " has a sentence-fragment summary: " + summary);
                end
            end
        end

        function transformDomainDocumentationStatesContractsAndDefaults(testCase)
            hydrostatic = testCase.generatedTransformPage("wvtransformhydrostatic","index.md");
            barotropic = testCase.generatedTransformPage("wvtransformbarotropicqg","index.md");
            expectations = {
                "x.md", ["Nx`-by-1","meters","endpoint"]
                "spatialmatrixsize.md", ["[Nx Ny Nz]","[Nx Ny]"]
                "spectralmatrixsize.md", ["[Nj Nkl]","[1 Nkl]"]
                "latitude.md", ["default is `33`","degrees north"]
                "planetaryradius.md", ["`6.371e6` m","beta"]
                "rotationrate.md", ["7.2921\times10^{-5}","\mathrm{rad\,s^{-1}}","inertialPeriod"]
                "g.md", ["9.81","\mathrm{m\,s^{-2}}"]
                "rho0.md", ["1025","\mathrm{kg\,m^{-3}}"]
                "shouldantialias.md", ["default is `true`","two-thirds"]
                "shouldusetruenomotionprofile.md", ["default is `false`","rhoFunction"]
                "h_0.md", ["Nj`-by-1","meters"]
                "h_pm.md", ["[Nj Nkl]","Omega"]
                "lr2.md", ["g h_0/f^2","square meters"]
                };
            for iExpectation = 1:size(expectations,1)
                page = testCase.generatedTransformPage("wvtransformhydrostatic",expectations{iExpectation,1});
                for phrase = expectations{iExpectation,2}
                    testCase.verifySubstring(page,phrase,expectations{iExpectation,1});
                end
            end
            testCase.verifySubstring(hydrostatic,"[`xyzGrid`]");
            testCase.verifySubstring(barotropic,"[`xyGrid`]");
            testCase.verifySubstring(testCase.generatedTransformPage("wvtransformbarotropicqg","h.md"),"default is `0.8` m");
            testCase.verifySubstring(testCase.generatedTransformPage("wvtransformconstantstratification","n0.md"),"5.2\times10^{-3}");
            testCase.verifySubstring(testCase.generatedTransformPage("wvtransformconstantstratification","n0.md"),"\mathrm{rad\,s^{-1}}");
        end

        function transformDocumentationExamplesExecute(testCase)
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5],N0=5.2e-3,latitude=45,shouldAntialias=true);
            wvt.initWithWaveModes(kMode=1,lMode=0,j=1,u=0.01,sign=1);

            [X,Y,Z] = wvt.xyzGrid;
            [K,L,J] = wvt.kljGrid;
            testCase.verifySize(X,wvt.spatialMatrixSize);
            testCase.verifySize(Y,wvt.spatialMatrixSize);
            testCase.verifySize(Z,wvt.spatialMatrixSize);
            testCase.verifySize(K,wvt.spectralMatrixSize);
            testCase.verifySize(L,wvt.spectralMatrixSize);
            testCase.verifySize(J,wvt.spectralMatrixSize);

            [u,v,eta] = wvt.variableWithName('u','v','eta');
            [uPoint,vPoint] = wvt.variableAtPositionWithName(wvt.x(2),wvt.y(2),wvt.z(2),'u','v',interpolationMethod='spline');
            testCase.verifySize(u,wvt.spatialMatrixSize);
            testCase.verifySize(v,wvt.spatialMatrixSize);
            testCase.verifySize(eta,wvt.spatialMatrixSize);
            testCase.verifyTrue(all(isfinite([uPoint vPoint])));
            testCase.verifyGreaterThanOrEqual(wvt.totalEnergyOfFlowComponent(wvt.waveComponent),0);

            N2 = @(z)(5.2e-3)^2*ones(size(z));
            hydrostatic = WVTransformHydrostatic([4000 3000 1000],[8 6 5],N2Function=N2,latitude=45,shouldAntialias=false);
            testCase.verifyEqual(hydrostatic.volumeIntegral(ones(hydrostatic.spatialMatrixSize)),hydrostatic.Lz,RelTol=1e-12);

            dudx = wvt.diffX(u);
            dudz = wvt.diffZF(u);
            antiderivative = wvt.intZG(dudz);
            testCase.verifySize(dudx,wvt.spatialMatrixSize);
            testCase.verifySize(dudz,wvt.spatialMatrixSize);
            testCase.verifySize(antiderivative,wvt.spatialMatrixSize);

            radial = wvt.transformToRadialWavenumber(wvt.Apm_TE_factor.*(abs(wvt.Ap).^2+abs(wvt.Am).^2));
            testCase.verifyEqual(size(radial,1),wvt.Nj);
            testCase.verifyTrue(all(isfinite(radial),"all"));

            resized = wvt.waveVortexTransformWithResolution([10 8 7]);
            explicit = wvt.waveVortexTransformWithExplicitAntialiasing;
            testCase.verifyClass(resized,"WVTransformConstantStratification");
            testCase.verifyClass(explicit,"WVTransformConstantStratification");
            testCase.verifyEqual(resized.spatialMatrixSize,[10 8 7]);
            testCase.verifyFalse(explicit.shouldAntialias);

            qg = WVTransformBarotropicQG([4000 3000],[8 6],h=0.8,latitude=45);
            [Xqg,Yqg] = qg.xyGrid;
            [Kqg,Lqg] = qg.klGrid;
            qg.setSSH(@(x,y)0.1*cos(2*pi*x/qg.Lx));
            A0 = qg.transformQGPVToWaveVortex(qg.qgpv);
            testCase.verifySize(Xqg,qg.spatialMatrixSize);
            testCase.verifySize(Yqg,qg.spatialMatrixSize);
            testCase.verifySize(Kqg,qg.spectralMatrixSize);
            testCase.verifySize(Lqg,qg.spectralMatrixSize);
            testCase.verifySize(A0,qg.spectralMatrixSize);
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

        function forcingReferenceUsesReviewedTopics(testCase)
            expectedOrder = [
                "Create the forcing"
                "Inspect forcing configuration"
                "Configure forcing"
                "Inspect forcing or damping scales"
                "Evaluate prescribed forcing"
                "Generate forcing inputs"
                ];
            indexes = [
                fullfile("wvforcing","index.md")
                fullfile("wvnonlinearadvection","index.md")
                fullfile("wvbottomfrictionlinear","index.md")
                fullfile("wvbottomfrictionquadratic","index.md")
                fullfile("wvfixedamplitudeforcing","index.md")
                fullfile("wvnarrowbandgeostrophicforcing","index.md")
                fullfile("wvbetaplanepvadvection","index.md")
                fullfile("wvpseudotopographicwavegeneration","index.md")
                fullfile("closures","wvadaptivedamping","index.md")
                fullfile("closures","wvverticaldiffusivity","index.md")
                fullfile("closures","wvhorizontaldamping","index.md")
                fullfile("closures","wvverticaldamping","index.md")
                fullfile("closures","wvthermaldamping","index.md")
                fullfile("closures","wvantialiasing","index.md")
                ];
            for relativePath = indexes'
                page = testCase.generatedForcingPage(relativePath);
                topics = testCase.topLevelTopics(page);
                testCase.verifyEqual(topics,expectedOrder(ismember(expectedOrder,topics)),relativePath);
                testCase.verifyFalse(contains(page,newline + "+ Other" + newline),relativePath);
                testCase.verifyNoEmptyTopicBranches(page,relativePath);
            end

            expectedUserMembers = {
                fullfile("wvforcing","index.md"), ["wvt" "name" "forcingType" "isClosure" "priority"]
                fullfile("wvbottomfrictionlinear","index.md"), ["r" "r_scaled"]
                fullfile("wvbottomfrictionquadratic","index.md"), ["Cd" "cd"]
                fullfile("wvfixedamplitudeforcing","index.md"), ["A0_indices" "A0bar" "Ap_indices" "Apbar" "Am_indices" "Ambar" "setWaveForcingCoefficients" "setGeostrophicForcingCoefficients" "setNarrowBandGeostrophicForcing"]
                fullfile("wvnarrowbandgeostrophicforcing","index.md"), ["r" "k_r" "k_f" "j_f" "u_rms" "initialPV" "modelSpectrum" "A0_indices" "A0bar"]
                fullfile("closures","wvadaptivedamping","index.md"), ["damp" "k_no_damp" "k_damp" "j_no_damp" "j_damp" "dampingTimeScale"]
                fullfile("closures","wvverticaldiffusivity","index.md"), ["kappa_z" "shouldForceMeanDensityAnomaly"]
                fullfile("closures","wvhorizontaldamping","index.md"), ["nu" "kappa"]
                fullfile("closures","wvverticaldamping","index.md"), ["nu" "kappa"]
                fullfile("closures","wvthermaldamping","index.md"), ["alpha" "alpha_scaled"]
                fullfile("closures","wvantialiasing","index.md"), ["Nj" "effectiveHorizontalGridResolution" "effectiveJMax"]
                fullfile("wvpseudotopographicwavegeneration","index.md"), ["topographicHeight" "barotropicVelocityAmplitude" "frequency" "darwinSymbol" "rampDuration" "startTime" "shouldAvoidAdaptiveDamping" "maximumForcedHorizontalWavenumber" "maximumForcedVerticalMode" "barotropicVelocityAtTime" "bottomVelocityAtTime" "spectralGenerationMask" "goffAbyssalHillTopography"]
                };
            for iExpectation = 1:size(expectedUserMembers,1)
                relativePath = expectedUserMembers{iExpectation,1};
                page = testCase.generatedForcingPage(relativePath);
                userTopics = extractBetween(page,"## Topics","## Developer Topics");
                for member = expectedUserMembers{iExpectation,2}
                    testCase.verifySubstring(userTopics,"[`" + member + "`]",relativePath + " hides " + member + ".");
                end
            end

            base = testCase.generatedForcingPage(fullfile("wvforcing","index.md"));
            baseUserTopics = extractBetween(base,"## Topics","## Developer Topics");
            baseDeveloperTopics = extractAfter(base,"## Developer Topics");
            evaluationHooks = [
                "addHydrostaticSpatialForcing"
                "addNonhydrostaticSpatialForcing"
                "addPotentialVorticitySpatialForcing"
                "addSpectralForcing"
                "setSpectralForcing"
                "setSpectralAmplitude"
                ];
            for hook = evaluationHooks'
                testCase.verifyFalse(contains(baseUserTopics,"[`" + hook + "`]"));
                testCase.verifySubstring(baseDeveloperTopics,"[`" + hook + "`]");
            end
            antialiasing = testCase.generatedForcingPage(fullfile("closures","wvantialiasing","index.md"));
            testCase.verifySubstring(extractAfter(antialiasing,"## Developer Topics"),"[`M`]");
        end

        function forcingBaseDocumentsStagesAndConstructor(testCase)
            constructor = testCase.generatedForcingPage(fullfile("wvforcing","wvforcing.md"));
            testCase.verifySubstring(constructor,"self = WVForcing(wvt,name,forcingType)");
            for parameter = ["wvt" "name" "forcingType"]
                testCase.verifySubstring(constructor,"`" + parameter + "`");
            end
            testCase.verifySubstring(constructor,"`self`");
            testCase.verifyFalse(contains(constructor,"WVNonlinearFluxOperation"));

            forcingType = testCase.generatedForcingPage(fullfile("wvforcing","forcingtype.md"));
            stageMappings = {
                "HydrostaticSpatial", "addHydrostaticSpatialForcing"
                "NonhydrostaticSpatial", "addNonhydrostaticSpatialForcing"
                "PVSpatial", "addPotentialVorticitySpatialForcing"
                "Spectral", "addSpectralForcing"
                "PVSpectral", "addPotentialVorticitySpectralForcing"
                "SpectralAmplitude", "setSpectralForcing"
                "PVSpectralAmplitude", "setPotentialVorticitySpectralForcing"
                };
            for iMapping = 1:size(stageMappings,1)
                testCase.verifySubstring(forcingType,"`" + stageMappings{iMapping,1} + "`");
                testCase.verifySubstring(forcingType,"`" + stageMappings{iMapping,2} + "`");
            end
            testCase.verifySubstring(forcingType,"restores the constrained coefficient values exactly");
            priority = testCase.generatedForcingPage(fullfile("wvforcing","priority.md"));
            testCase.verifySubstring(priority,"from 0 first to 255 last");
            testCase.verifySubstring(priority,"same evaluation stage");
        end

        function forcingLandingPagesCoverGeneratedClasses(testCase)
            forcingLanding = string(fileread(fullfile(testCase.repositoryRoot,"Documentation", ...
                "WebsiteDocumentation","classes","forcing","index.md")));
            closureLanding = string(fileread(fullfile(testCase.repositoryRoot,"Documentation", ...
                "WebsiteDocumentation","classes","forcing","closures","index.md")));
            forcingClasses = ["WVNonlinearAdvection" "WVBottomFrictionLinear" ...
                "WVBottomFrictionQuadratic" "WVFixedAmplitudeForcing" ...
                "WVNarrowBandGeostrophicForcing" "WVBetaPlanePVAdvection" ...
                "WVPseudoTopographicWaveGeneration"];
            closureClasses = ["WVAdaptiveDamping" "WVHorizontalDamping" ...
                "WVVerticalDamping" "WVVerticalDiffusivity" ...
                "WVThermalDamping" "WVAntialiasing"];
            for className = forcingClasses
                testCase.verifyEqual(count(forcingLanding,"[`" + className + "`]"),1,className);
            end
            for className = closureClasses
                testCase.verifyEqual(count(closureLanding,"[`" + className + "`]"),1,className);
            end
            testCase.verifySubstring(forcingLanding,"| Class | Purpose |");
            testCase.verifySubstring(closureLanding,"| Class | Purpose |");
            testCase.verifyFalse(contains(forcingLanding,"Supported transforms"));
            testCase.verifyFalse(contains(forcingLanding,"Principal controls"));
            testCase.verifyFalse(contains(closureLanding,"Supported transforms"));
            testCase.verifyFalse(contains(closureLanding,"Principal controls"));
            testCase.verifySubstring(closureLanding,"Transform-level antialiasing remains the efficient default");
            sidecarFolder = fullfile(testCase.repositoryRoot,"Forcing","detailedDescriptions");
            remainingSidecars = dir(fullfile(sidecarFolder,"*.md"));
            testCase.verifyEmpty(remainingSidecars);
        end

        function forcingDetailedDescriptionsRemainOnClassPages(testCase)
            expectations = {
                "wvbottomfrictionlinear", [
                    "bottom quadrature weight"
                    "r_\mathrm{scaled}"
                    "Comparing this with quadratic drag"
                    "For both nonhydrostatic and hydrostatic transforms"
                    "and for quasigeostrophic transforms"
                    "### Example"
                    "wvt.addForcing(WVBottomFrictionLinear(wvt,r=1/(200*86400)))"
                    ]
                "wvbottomfrictionquadratic", [
                    "quadrature weight"
                    "4000 m reference depth"
                    "Comparing quadratic and linear drag"
                    "hydrostatic and nonhydrostatic transforms"
                    "and for quasigeostrophic transforms"
                    "### Example"
                    "wvt.addForcing(WVBottomFrictionQuadratic(wvt,Cd=0.001))"
                    ]
                "wvfixedamplitudeforcing", [
                    "keeps selected wave-vortex coefficients at prescribed amplitudes"
                    "participate in nonlinear interactions"
                    "participate in all the nonlinear dynamics"
                    "spectral-amplitude forcing"
                    "prescribed coefficient values"
                    "removes a degree of freedom"
                    "closure's damping range"
                    "### Example"
                    ]
                "wvnonlinearadvection", [
                    "nonhydrostatic transforms"
                    "hydrostatic transforms"
                    "quasigeostrophic transforms"
                    "\mathcal{S}_w"
                    "\mathcal{S}_\eta"
                    "\mathcal{S}_\mathrm{qgpv}"
                    "installs this forcing by default"
                    ]
                };
            for iExpectation = 1:size(expectations,1)
                classFolder = expectations{iExpectation,1};
                overview = testCase.generatedForcingPage(fullfile(classFolder,"index.md"));
                for expectedText = expectations{iExpectation,2}'
                    testCase.verifySubstring(overview,expectedText,classFolder + " lost authored detail.");
                end

                constructor = testCase.generatedForcingPage(fullfile(classFolder,classFolder + ".md"));
                expectedLink = "](/classes/forcing/" + classFolder + "/)";
                testCase.verifySubstring(constructor,expectedLink,classFolder + " constructor does not link to its overview.");
            end
            fixedAmplitude = testCase.generatedForcingPage(fullfile("wvfixedamplitudeforcing","index.md"));
            testCase.verifyFalse(contains(fixedAmplitude,"recording the tendency that must be cancelled"));
        end

        function coefficientClaimsMatchTransformState(testCase)
            wvt = WVTransformConstantStratification([40e3 30e3 2e3],[8 6 5],N0=5.2e-3,latitude=45);
            testCase.verifyEqual(size(wvt.Ap),size(wvt.Am));
            testCase.verifyEqual(size(wvt.Ap),size(wvt.A0));
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("Ap").units),"m s-1");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("Am").units),"m s-1");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("A0").units),"m2 s-1");
            testCase.verifyEqual(string(wvt.propertyAnnotationWithName("A0t").units),"m2 s-1");

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

        function page = generatedForcingPage(testCase,relativePath)
            page = string(fileread(fullfile(testCase.repositoryRoot,"docs","classes","forcing",relativePath)));
        end

        function topics = topLevelTopics(~,page)
            topicBlock = extractBetween(page,"## Topics","## Developer Topics");
            tokens = regexp(topicBlock,'(?m)^\+ (?<name>[^\r\n]+)$','names');
            topics = string({tokens.name})';
        end

        function summary = memberSummary(testCase,page,name)
            lines = splitlines(page);
            candidate = lines(startsWith(strtrim(lines),"+ [`" + name + "`]"));
            testCase.assertNumElements(candidate,1,"Expected exactly one topic-list entry for " + name + ".");
            summary = strtrim(regexprep(candidate,'^.*\)\s*',''));
        end

        function verifyTextOrder(testCase,text,tokens)
            positions = zeros(size(tokens));
            for iToken = 1:numel(tokens)
                location = strfind(text,tokens(iToken));
                testCase.assertNotEmpty(location,"Missing " + tokens(iToken));
                positions(iToken) = location(1);
            end
            testCase.verifyGreaterThan(diff(positions),zeros(1,numel(tokens)-1), ...
                "Items are not in the requested order: " + strjoin(tokens," → "));
        end

        function verifyNoEmptyTopicBranches(testCase,page,label)
            lines = splitlines(extractBetween(page,"## Topics","## Developer Topics"));
            for iLine = 1:numel(lines)
                topic = regexp(lines(iLine),'^(?<indent> *)\+ (?<body>.+)$','names','once');
                if isempty(topic) || contains(string(topic.body),"[`")
                    continue
                end
                depth = strlength(string(topic.indent));
                hasLinkedMember = false;
                for iChild = (iLine+1):numel(lines)
                    child = regexp(lines(iChild),'^(?<indent> *)\+ (?<body>.+)$','names','once');
                    if isempty(child)
                        continue
                    end
                    if strlength(string(child.indent)) <= depth
                        break
                    end
                    hasLinkedMember = hasLinkedMember || contains(string(child.body),"[`");
                end
                testCase.verifyTrue(hasLinkedMember, ...
                    label + " has an empty topic branch: " + string(topic.body));
            end
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
