classdef TestPortableVariableCatalog < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function testGeneratedCatalogMatchesCommittedFiles(testCase)
            repositoryRoot = TestPortableVariableCatalog.repositoryRoot();
            originalPath = path;
            pathCleanup = onCleanup(@() path(originalPath));
            addpath(fullfile(repositoryRoot,"tools"));
            temporaryRoot = string(tempname);
            cleanup = onCleanup(@() rmdir(temporaryRoot,"s"));
            supplementFolder = fullfile(temporaryRoot,"PortableRuntime","contracts");
            mkdir(supplementFolder);
            copyfile(fullfile(repositoryRoot,"PortableRuntime","contracts", ...
                "portable-variable-supplement-v1.json"),supplementFolder);

            originalDirectory = string(pwd);
            directoryCleanup = onCleanup(@() cd(originalDirectory));
            cd(tempdir);
            result = generatePortableVariableCatalog(repositoryRoot=temporaryRoot);
            clear directoryCleanup

            testCase.verifyEqual(fileread(result.catalogPath),fileread(fullfile( ...
                repositoryRoot,"PortableRuntime","contracts", ...
                "portable-variable-catalog-v1.json")));
            testCase.verifyEqual(fileread(result.headerPath),fileread(fullfile( ...
                repositoryRoot,"PortableRuntime","include", ...
                "WaveVortexRuntime","generated", ...
                "WVPortableVariableCatalog.hpp")));
            testCase.verifyEqual(fileread(result.contractHeaderPath),fileread(fullfile( ...
                repositoryRoot,"PortableRuntime","include","WaveVortexRuntime", ...
                "generated","WVPortableVariableContracts.hpp")));
            testCase.verifyEqual(fileread(result.documentationPath),fileread(fullfile( ...
                repositoryRoot,"PortableRuntime","VARIABLES.md")));
            clear pathCleanup
            clear cleanup
        end

        function testCatalogIdentityAndAnnotationMetadata(testCase)
            repositoryRoot = TestPortableVariableCatalog.repositoryRoot();
            catalog = jsondecode(fileread(fullfile(repositoryRoot, ...
                "PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            variables = catalog.variables;
            testCase.verifyEqual(string(catalog.schema),"portable-variable-catalog-v1");
            testCase.verifyEqual([variables.ordinal],0:75);
            testCase.verifyEqual(numel(unique(string({variables.name}))),76);

            u = variables(string({variables.name}) == "u");
            testCase.verifyEqual(string(u.dimensions),["x";"y";"z"]);
            testCase.verifyEqual(string(u.units),"m s-1");
            testCase.verifyEqual(string(u.description), ...
                "x-component of the fluid velocity");
            testCase.verifyFalse(u.isComplex);
            testCase.verifyEqual(double(u.primitiveDependencyMask),1);
            testCase.verifyEqual(string(u.netCDFAttributes.name),"standard_name");
            testCase.verifyEqual(string(u.netCDFAttributes.value), ...
                "eastward_sea_water_velocity");

            Ap = variables(string({variables.name}) == "Ap");
            testCase.verifyEqual(string(Ap.dimensions),["j";"kl"]);
            testCase.verifyTrue(Ap.isComplex);
            testCase.verifyFalse(Ap.isVariableWithLinearTimeStep);
            testCase.verifyTrue(Ap.isVariableWithNonlinearTimeStep);
        end

        function testLegacyMetadataAndGraphFailures(testCase)
            root = TestPortableVariableCatalog.repositoryRoot();
            catalog = jsondecode(fileread(fullfile(root,"PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            legacy = jsondecode(fileread(fullfile(root,"UnitTests","fixtures","portable-variable-legacy-v1.json")));
            testCase.verifyEqual(catalog.variables(1:23),legacy.variables);
            validatePortableVariableCatalog(catalog);
            testCase.verifyNumElements(catalog.configurations,12);
            testCase.verifyNumElements(catalog.contracts,666);
            broken = catalog;
            broken.contracts(1).dependencies = {broken.contracts(1).metadata.name};
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableDependency");
            broken = catalog;
            broken.contracts(1).dependencies = {'customUnknown'};
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableDependency");
            broken = catalog;
            broken.contracts(2) = broken.contracts(1);
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
            broken = catalog;
            broken.contracts(1).metadata.units = '';
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
            broken = catalog;
            broken.contracts(1).intermediateMask = 2^30;
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
            broken = catalog;
            broken.exclusions(2) = broken.exclusions(1);
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
            broken = catalog;
            broken.exclusions(1).name = 'u';
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
            broken = catalog;
            broken.contracts(1).metadata.samplingModes = {'positions'};
            testCase.verifyError(@()validatePortableVariableCatalog(broken),"WaveVortexModel:InvalidPortableVariableCatalog");
        end

        function testConfigurationAndConditionalDependencies(testCase)
            root = TestPortableVariableCatalog.repositoryRoot();
            catalog = jsondecode(fileread(fullfile(root,"PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            rows = catalog.contracts;
            names = arrayfun(@(r)string(r.metadata.name),rows);
            barotropic = rows(names=="u" & string({rows.configuration}).'=="barotropic-aa0");
            testCase.verifyEqual(string(barotropic.metadata.dimensions),["x";"y"]);
            hydrostatic = rows(names=="eta_true" & string({rows.configuration}).'=="hydrostatic-aa0");
            testCase.verifyEqual(string(hydrostatic.trueProfileDependencies),"rho_nm");
            testCase.verifyEqual(string(hydrostatic.dependencies),"rho_total");
            ape = rows(names=="ape" & string({rows.configuration}).'=="hydrostatic-aa0");
            testCase.verifyEqual(string(ape.dependencies),"eta_true");
            testCase.verifyEqual(string(ape.trueProfileDependencies),"rho_nm");
            testCase.verifyEqual(string(catalog.policy.defaultNoMotionSolver),"dampedLeastSquares");
            testCase.verifyTrue(catalog.policy.shouldUseTrueNoMotionProfileDefault);
            testCase.verifyEqual(string(catalog.policy.noMotionSolvers),["dampedLeastSquares";"lsqnonlin";"fminsearch"]);
            testCase.verifyEqual(string(hydrostatic.runtimeStatus),"implemented");
            forcing = rows(names=="Fqgpv_portable_catalog_forcing" & string({rows.configuration}).'=="barotropic-aa0");
            testCase.verifyEqual(string(forcing.metadata.units),"s-2");
            testCase.verifyEqual(string(forcing.authority),"forcing-instance-template");
            testCase.verifyEqual(string(forcing.runtimeStatus),"implemented");
            testCase.verifyTrue(ismember("totalEnstrophy",string({catalog.exclusions.name})));
        end

        function testEvaluationLoopsDoNotCompareFieldNames(testCase)
            repositoryRoot = TestPortableVariableCatalog.repositoryRoot();
            source = fileread(fullfile(repositoryRoot,"PortableRuntime","src", ...
                "WVFieldEvaluationService.cpp"));
            testCase.verifyEmpty(regexp(source, ...
                'request\.fieldName\s*(==|!=)',"once"));
            testCase.verifyNotEmpty(strfind(source, ...
                "findExecutablePortableVariable(request.fieldName)"));
        end

        function testHydrostaticAndNonhydrostaticAnnotationsAgree(testCase)
            fieldNames = {'u','v','w','eta','pi','p','psi','qgpv', ...
                'rho_e','rho_total','rho_bar','zeta_x','zeta_y','zeta_z', ...
                'ssu','ssv','ssh','energy','uvMax','wMax'};
            hydrostatic = WVTransformHydrostatic.classDefinedOperationForKnownVariable( ...
                fieldNames{:},spatialDimensionNames={'x','y','z'});
            nonhydrostatic = WVTransformBoussinesq.classDefinedOperationForKnownVariable( ...
                fieldNames{:},spatialDimensionNames={'x','y','z'});
            for iField = 1:numel(fieldNames)
                first = hydrostatic(iField).outputVariables;
                second = nonhydrostatic(iField).outputVariables;
                testCase.verifyEqual(first.name,second.name);
                testCase.verifyEqual(first.dimensions,second.dimensions);
                testCase.verifyEqual(first.units,second.units);
                testCase.verifyEqual(first.description,second.description);
                testCase.verifyEqual(first.isComplex,second.isComplex);
                testCase.verifyEqual(first.isVariableWithLinearTimeStep, ...
                    second.isVariableWithLinearTimeStep);
                testCase.verifyEqual(first.isVariableWithNonlinearTimeStep, ...
                    second.isVariableWithNonlinearTimeStep);
            end
        end
    end

    methods (Static, Access=private)
        function root = repositoryRoot()
            root = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end
end
