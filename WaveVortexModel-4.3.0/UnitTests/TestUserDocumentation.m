classdef TestUserDocumentation < matlab.unittest.TestCase
    properties
        repositoryRoot (1,1) string
        canonicalRoot (1,1) string
    end

    methods (TestMethodSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.canonicalRoot = fullfile(testCase.repositoryRoot,"Documentation","WebsiteDocumentation");
        end
    end

    methods (Test,TestTags="full")
        function readmeAndHomepageShareExecutableQuickStart(testCase)
            readmeBlock = testCase.quickStartBlock(fullfile(testCase.repositoryRoot,"README.md"));
            homepageBlock = testCase.quickStartBlock(fullfile(testCase.canonicalRoot,"index.md"));
            testCase.verifyEqual(readmeBlock,homepageBlock);

            eval(char(readmeBlock));
            testCase.verifyClass(wvt,"WVTransformConstantStratification");
            testCase.verifyClass(model,"WVModel");
            testCase.verifyFalse(model.isDynamicsLinear);
            testCase.verifyClass(wvt.forcingWithName("adaptive damping"),"WVAdaptiveDamping");
            testCase.verifyEqual(model.t,600);
            testCase.verifyTrue(any(abs(u)>0,"all"));
            testCase.verifySize(v,wvt.spatialMatrixSize);
            testCase.verifySize(w,wvt.spatialMatrixSize);
            testCase.verifyTrue(all(isfinite([omega k l])));
            testCase.verifyFalse(contains(readmeBlock,"shouldUseLinearDynamics"));
            testCase.verifyFalse(contains(readmeBlock,"shouldShowIntegrationDiagnostics"));
            testCase.verifyFalse(contains(readmeBlock,"..."));
        end

        function readmePointsToHostedDocumentation(testCase)
            readme = testCase.readFile(fullfile(testCase.repositoryRoot,"README.md"));
            linkLocations = strfind(readme,"https://wavevortexmodel.org");
            testCase.verifyLessThan(linkLocations(1),500, ...
                "The canonical documentation link must be prominent.");
            requiredRoutes = [
                "https://wavevortexmodel.org/installation"
                "https://wavevortexmodel.org/users-guide/"
                "https://wavevortexmodel.org/users-guide/supported-features.html"
                "https://wavevortexmodel.org/classes/transforms/"
                "https://wavevortexmodel.org/classes/wvmodel/"
                "https://wavevortexmodel.org/mathematical-introduction/"
                ];
            for route = requiredRoutes'
                testCase.verifySubstring(readme,route);
            end
            for citation = ["10.5281/zenodo.4037401" "10.1017/jfm.2020.995" "10.48550/arXiv.2403.20269"]
                testCase.verifySubstring(readme,citation);
            end
        end

        function onboardingUsesCurrentInterfaces(testCase)
            homepage = testCase.readCanonical("index.md");
            introduction = testCase.readCanonical(fullfile("users-guide","introduction.md"));
            transformGuide = testCase.readCanonical(fullfile("users-guide","using-the-wvtransform.md"));
            persistenceGuide = testCase.readCanonical(fullfile("users-guide","reading-and-writing-to-file.md"));
            combined = homepage + introduction + transformGuide + persistenceGuide;

            for className = [
                    "WVTransformConstantStratification"
                    "WVTransformHydrostatic"
                    "WVTransformBoussinesq"
                    "WVTransformStratifiedQG"
                    "WVTransformBarotropicQG"
                    "WVModel"
                    ]'
                testCase.verifySubstring(combined,className);
            end
            for entryPoint = [
                    "initWithWaveModes"
                    "initWithUVEta"
                    "variableWithName"
                    "variableAtPositionWithName"
                    "waveVortexTransformFromFile"
                    "modelFromFile"
                    ]'
                testCase.verifySubstring(combined,entryPoint);
            end
            testCase.verifySubstring(transformGuide,"`linear` or `spline`");
            testCase.verifySubstring(transformGuide,"isHydrostatic=true");

            nonlinearStatement = "By default, it integrates nonlinear interactions among the resolved flow components.";
            linearStatement = "Analytical linear evolution is also available when those nonlinear interactions should be omitted.";
            nonlinearLocation = strfind(introduction,nonlinearStatement);
            linearLocation = strfind(introduction,linearStatement);
            testCase.verifyNumElements(nonlinearLocation,1);
            testCase.verifyNumElements(linearLocation,1);
            testCase.verifyLessThan(nonlinearLocation,linearLocation);
            testCase.verifyFalse(contains(persistenceGuide,"shouldUseLinearDynamics"));
        end

        function handAuthoredExamplesAvoidUnnecessaryContinuations(testCase)
            markdownFiles = dir(fullfile(testCase.canonicalRoot,"**","*.md"));
            paths = [fullfile(testCase.repositoryRoot,"README.md"); string(fullfile({markdownFiles.folder},{markdownFiles.name})).'];
            for path = paths'
                page = testCase.readFile(path);
                blocks = regexp(page,'(?s)```matlab\r?\n(?<code>.*?)```','names');
                for iBlock = 1:numel(blocks)
                    testCase.verifyEmpty(regexp(string(blocks(iBlock).code),'(?m)\.\.\.\s*$','once'), ...
                        "Unnecessary MATLAB continuation remains in " + path);
                end
            end
        end

        function canonicalGuidanceAvoidsObsoleteProse(testCase)
            markdownFiles = dir(fullfile(testCase.canonicalRoot,"**","*.md"));
            contents = strings(numel(markdownFiles)+1,1);
            contents(1) = testCase.readFile(fullfile(testCase.repositoryRoot,"README.md"));
            for iFile = 1:numel(markdownFiles)
                contents(iFile+1) = testCase.readFile(fullfile(markdownFiles(iFile).folder,markdownFiles(iFile).name));
            end
            allGuidance = join(contents,newline);
            forbidden = [
                "InternalWaveModel"
                "exact interpolation"
                "external-wave interpolation"
                "| Status |"
                "Status definitions"
                "Stable"
                "Experimental"
                "Deprecated"
                "v4.2.1"
                ];
            for phrase = forbidden'
                testCase.verifyFalse(contains(allGuidance,phrase), ...
                    "Obsolete or project-oriented phrase remains: " + phrase);
            end
        end

        function netCDFConventionsAreSubstantive(testCase)
            page = testCase.readCanonical(fullfile("users-guide","netcdf-output.md"));
            for heading = ["## Scientific variables" "## CF conventions" "## Discovery metadata" "## Related guidance"]
                testCase.verifySubstring(page,heading);
            end
            testCase.verifySubstring(page,"WVVariableAnnotation");
            testCase.verifySubstring(page,"does not claim");
            testCase.verifySubstring(page,"responsibility of the file producer");
            testCase.verifyFalse(contains(page,"## Global attributes"+newline+newline+"## Dimensions"));
        end

        function compiledExecutionNavigationAndRoutesAreStable(testCase)
            overview = testCase.readCanonical(fullfile("compiled-execution","index.md"));
            matlabPreview = testCase.readCanonical(fullfile( ...
                "compiled-execution","compiled-matlab-backend.md"));
            standalone = testCase.readCanonical(fullfile( ...
                "compiled-execution","standalone-portable-runtime.md"));
            userGuide = testCase.readCanonical(fullfile("users-guide","index.md"));
            benchmarks = testCase.readCanonical("benchmarks.md");

            testCase.verifyMatches(overview,'(?m)^title: Compiled execution$')
            testCase.verifyMatches(overview,'(?m)^nav_order: 9\.5$')
            testCase.verifyMatches(overview,'(?m)^permalink: /compiled-execution$')
            testCase.verifyMatches(matlabPreview, ...
                '(?m)^permalink: /users-guide/compiled-preview\.html$')
            testCase.verifyMatches(standalone, ...
                '(?m)^permalink: /users-guide/portable-runtime\.html$')
            testCase.verifyMatches(matlabPreview,'(?m)^parent: Compiled execution$')
            testCase.verifyMatches(standalone,'(?m)^parent: Compiled execution$')
            testCase.verifySubstring(userGuide,"[Compiled execution](/compiled-execution)")
            testCase.verifySubstring(standalone, ...
                'WVModel.writePortableRunRequest("run.json","initial-condition.nc",finalTime=86400);')
            testCase.verifySubstring(standalone,"wave-vortex-run --request run.json")
            testCase.verifySubstring(overview,"Compiled MATLAB backend preview")
            testCase.verifySubstring(overview,"Standalone portable runtime")
            testCase.verifySubstring(benchmarks,"frozen to the accepted issue #312")
            testCase.verifySubstring(benchmarks,"integration-only runtime and total peak process memory")
            testCase.verifyFalse(isfile(fullfile(testCase.canonicalRoot, ...
                "users-guide","compiled-preview.md")))
            testCase.verifyFalse(isfile(fullfile(testCase.canonicalRoot, ...
                "users-guide","portable-runtime.md")))
        end

        function websiteUsesModernTypographyAndSearch(testCase)
            expectedFont = 'font-family: "Avenir Next", Avenir, "Helvetica Neue", Helvetica, Arial, sans-serif;';
            for root = [testCase.canonicalRoot fullfile(testCase.repositoryRoot,"docs")]
                stylesheet = testCase.readFile(fullfile(root,"_sass","custom","custom.scss"));
                config = testCase.readFile(fullfile(root,"_config.yml"));
                testCase.verifySubstring(stylesheet,expectedFont);
                testCase.verifySubstring(stylesheet,"line-height: 1.65;");
                testCase.verifySubstring(stylesheet,"letter-spacing: -0.015em;");
                testCase.verifySubstring(stylesheet,".search-input");
                testCase.verifyMatches(config,"(?m)^search_enabled: true$");
            end
        end

        function canonicalDomainIsPreserved(testCase)
            canonicalCNAME = strtrim(testCase.readFile(fullfile(testCase.canonicalRoot,"CNAME")));
            generatedCNAME = strtrim(testCase.readFile(fullfile(testCase.repositoryRoot,"docs","CNAME")));
            testCase.verifyEqual(canonicalCNAME,"wavevortexmodel.org");
            testCase.verifyEqual(generatedCNAME,canonicalCNAME);
        end

        function renderingSensitiveMarkupIsWellFormed(testCase)
            transformGuide = testCase.readCanonical(fullfile("users-guide","using-the-wvtransform.md"));
            testCase.verifySubstring(transformGuide,'$$5 \leq \lvert\mathrm{latitude}\rvert \leq 85$$');
            testCase.verifyFalse(contains(transformGuide,'|\mathrm{latitude}|'));

            forcingGuide = testCase.readCanonical(fullfile("users-guide","adding-forcing.md"));
            equationStart = strfind(forcingGuide,"$$" + newline + "    \begin{align}");
            equationEnd = strfind(forcingGuide,"    \end{align}" + newline + "$$");
            explanatoryText = strfind(forcingGuide,"and implemented by overriding `addSpectralForcing`");
            testCase.verifyNumElements(equationStart,1);
            testCase.verifyNumElements(equationEnd,1);
            testCase.verifyNumElements(explanatoryText,1);
            testCase.verifyLessThan(equationStart,equationEnd);
            testCase.verifyLessThan(equationEnd,explanatoryText);
            testCase.verifySubstring(forcingGuide,'$$\psi = p/(\rho_0 f)$$');
            testCase.verifySubstring(forcingGuide,'$$\rho=- (\rho_0 f/g) \partial_z \psi$$');
            testCase.verifyEmpty(regexp(forcingGuide,'(?<![\\$])\$(?!\$)','once'));

            advancedGuide = testCase.readCanonical(fullfile("users-guide","reading-and-writing-to-file-advanced.md"));
            testCase.verifyMatches(advancedGuide,'(?m)^title: "Reading and writing files: advanced topics"$');
            testCase.verifyMatches(advancedGuide,'(?m)^parent: User guide$');

            landingPages = [
                "acknowledgements.md"
                "addons.md"
                fullfile("classes","index.md")
                fullfile("compiled-execution","index.md")
                fullfile("developers-guide","index.md")
                fullfile("mathematical-introduction","index.md")
                ];
            for landingPage = landingPages'
                testCase.verifyMatches(testCase.readCanonical(landingPage),'(?m)^# \S');
            end

            operationsGuide = testCase.readCanonical(fullfile("developers-guide","operations-and-variables.md"));
            propertyGuide = testCase.readCanonical(fullfile("developers-guide","property-and-method-types.md"));
            styleGuide = testCase.readCanonical(fullfile("developers-guide","style-guide.md"));
            changelog = testCase.readFile(fullfile(testCase.repositoryRoot,"CHANGELOG.md"));
            testCase.verifySubstring(operationsGuide,"`operationVariableNameMap`");
            testCase.verifySubstring(propertyGuide,"`WVTransform`");
            testCase.verifySubstring(styleGuide,"`WaveVortexTransformConstantStratification(Lxyz, Nxyz, N0, options)`");
            testCase.verifySubstring(changelog,"`WVForcing`");
            testCase.verifySubstring(changelog,"`WVObservingSystems`");
        end
    end

    methods (Access=private)
        function block = quickStartBlock(testCase,path)
            page = testCase.readFile(path);
            section = extractAfter(page,"## Quick start");
            block = extractBetween(section,"```matlab","```");
            testCase.assertNotEmpty(block,"No MATLAB quick-start block found in " + path);
            block = strtrim(block(1));
        end

        function contents = readCanonical(testCase,relativePath)
            contents = testCase.readFile(fullfile(testCase.canonicalRoot,relativePath));
        end

        function contents = readFile(~,path)
            contents = string(fileread(path));
        end
    end
end
