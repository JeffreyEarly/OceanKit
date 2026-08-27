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
            clear pathCleanup
            clear cleanup
        end

        function testCatalogIdentityAndAnnotationMetadata(testCase)
            repositoryRoot = TestPortableVariableCatalog.repositoryRoot();
            catalog = jsondecode(fileread(fullfile(repositoryRoot, ...
                "PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            variables = catalog.variables;
            testCase.verifyEqual(string(catalog.schema),"portable-variable-catalog-v1");
            testCase.verifyEqual([variables.ordinal],0:22);
            testCase.verifyEqual(numel(unique(string({variables.name}))),23);

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

        function testEvaluationLoopsDoNotCompareFieldNames(testCase)
            repositoryRoot = TestPortableVariableCatalog.repositoryRoot();
            source = fileread(fullfile(repositoryRoot,"PortableRuntime","src", ...
                "WVFieldEvaluationService.cpp"));
            testCase.verifyEmpty(regexp(source, ...
                'request\.fieldName\s*(==|!=)',"once"));
            testCase.verifyNotEmpty(strfind(source, ...
                "findPortableVariable(request.fieldName)")); %#ok<STRIFCND>
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
