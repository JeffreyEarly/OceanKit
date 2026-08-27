classdef TestPortableRunRequestWriter < matlab.unittest.TestCase
    properties (SetAccess=private)
        RepositoryRoot (1,1) string
        TemporaryFolder (1,1) string
        ModelPath (1,1) string
        BarotropicQGPath (1,1) string
    end

    methods (TestMethodSetup)
        function createFixture(testCase)
            testCase.RepositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryFolder = string(fixture.Folder);
            testCase.ModelPath = fullfile(testCase.TemporaryFolder,"model bundle Ω.nc");
            copyfile(testCase.fixturePath(),testCase.ModelPath);
            ncwriteatt(testCase.ModelPath,"/","portableFileIdentifier","primary");
            testCase.BarotropicQGPath = fullfile(testCase.TemporaryFolder,"qg model Ω.nc");
            testCase.createBarotropicQGFixture(testCase.BarotropicQGPath,[6 5],1,true);
            ncwriteatt(testCase.BarotropicQGPath,"/","portableFileIdentifier","qg-primary");
        end
    end

    methods (Test, TestTags="full")
        function writesMinimalV2RequestWithRuntimeDefaults(testCase)
            request = fullfile(testCase.TemporaryFolder,"minimal request.json");
            WVModel.writePortableRunRequest(request,testCase.ModelPath,finalTime=86400);

            document = jsondecode(fileread(request));
            testCase.verifyEqual(string(document.schemaIdentifier),"wave-vortex-run-request-v2")
            testCase.verifyEqual(document.schemaVersion,2)
            testCase.verifyEqual(string(fieldnames(document.integration)),"finalTime")
            testCase.verifyEqual(document.integration.finalTime,86400)
            testCase.verifyEmpty(fieldnames(document.execution))
            testCase.verifyEqual(string(document.output.policy),"append")

            firstBytes = testCase.fileBytes(request);
            WVModel.writePortableRunRequest(request,testCase.ModelPath,finalTime=86400);
            testCase.verifyEqual(testCase.fileBytes(request),firstBytes)
        end

        function writesOnlyExplicitV2Overrides(testCase)
            request = fullfile(testCase.TemporaryFolder,"overrides.json");
            WVModel.writePortableRunRequest(request,testCase.ModelPath,finalTime=86400, ...
                method="adaptive-rk45",initialStep=20,maximumStep=300, ...
                relativeTolerance=1e-7,absoluteToleranceScale=1e-9, ...
                fftProvider="reference",threads=1);

            document = jsondecode(fileread(request));
            testCase.verifyEqual(string(document.integration.method),"adaptive-rk45")
            testCase.verifyEqual(document.integration.initialStep,20)
            testCase.verifyEqual(document.integration.maximumStep,300)
            testCase.verifyEqual(document.integration.relativeTolerance,1e-7)
            testCase.verifyEqual(document.integration.absoluteToleranceScale,1e-9)
            testCase.verifyEqual(string(document.execution.fftProvider),"reference")
            testCase.verifyEqual(document.execution.threads,1)
        end

        function writesEveryAcceptedSchemaFormDeterministically(testCase)
            requestFolder = fullfile(testCase.TemporaryFolder,"requests with Ω");
            mkdir(requestFolder)
            copyfile(testCase.ModelPath,fullfile(requestFolder,"model bundle Ω.nc"));
            modelReference = "model bundle Ω.nc";

            requests = [
                testCase.writeRequest(requestFolder,"v1-fixed",1,"fixed-rk4",initialStep=.25)
                testCase.writeRequest(requestFolder,"v1-rk23",1,"adaptive-rk23",initialStep=.25,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6)
                testCase.writeRequest(requestFolder,"v2-fixed",2,"fixed-rk4",initialStep=.25)
                testCase.writeRequest(requestFolder,"v2-cfl-advective",2,"fixed-rk4",cfl=.25,timeStepConstraint="advective")
                testCase.writeRequest(requestFolder,"v2-cfl-oscillatory",2,"fixed-rk4",cfl=.25,timeStepConstraint="oscillatory")
                testCase.writeRequest(requestFolder,"v2-cfl-min",2,"fixed-rk4",cfl=.25,timeStepConstraint="min")
                testCase.writeRequest(requestFolder,"v2-rk23",2,"adaptive-rk23",initialStep=.25,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6)
                testCase.writeRequest(requestFolder,"v2-rk45",2,"adaptive-rk45",initialStep=.25,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6)
                testCase.writeRequest(requestFolder,"v2-rk78",2,"adaptive-rk78",initialStep=.25,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6)
                ];

            for iRequest = 1:numel(requests)
                document = jsondecode(fileread(requests(iRequest)));
                testCase.verifyEqual(string(document.modelFiles),modelReference)
                testCase.verifyEqual(string(fieldnames(document)), ...
                    ["schemaIdentifier";"schemaVersion";"modelFiles";"integration";"output";"execution";"report"])
                text = string(fileread(requests(iRequest)));
                for forbidden = ["transform","forcing","observer","schedule","restart","state"]
                    quote = string(char(34));
                    testCase.verifyFalse(contains(lower(text),quote+forbidden+quote))
                end
            end
        end

        function pathsMapsAndOrderFollowDocumentContract(testCase)
            folderA = fullfile(testCase.TemporaryFolder,"A Ω");
            folderB = fullfile(testCase.TemporaryFolder,"B Ω");
            mkdir(folderA)
            mkdir(folderB)
            for folder = [folderA folderB]
                copyfile(testCase.ModelPath,fullfile(folder,"first model.nc"));
                second = fullfile(folder,"second model.nc");
                copyfile(testCase.ModelPath,second);
                ncwriteatt(second,"/","portableFileIdentifier","secondary");
            end
            destinations = configureDictionary("string","string");
            destinations(["secondary","primary"]) = ["second output Ω.nc","first output.nc"];
            requestA = fullfile(folderA,"request.json");
            requestB = fullfile(folderB,"request.json");
            originalFolder = string(pwd);
            folderCleanup = onCleanup(@()cd(originalFolder));
            cd(tempdir)
            for request = [requestA requestB]
                WVModel.writePortableRunRequest(request,["second model.nc","first model.nc"], ...
                    schemaVersion=2,method="fixed-rk4",finalTime=4,initialStep=.25, ...
                    outputPolicy="create",destinations=destinations,fftProvider="reference",threads=1, ...
                    reportPath="report Ω.json");
            end
            clear folderCleanup

            testCase.verifyEqual(testCase.fileBytes(requestA),testCase.fileBytes(requestB))
            document = jsondecode(fileread(requestA));
            testCase.verifyEqual(string(document.modelFiles),["second model.nc";"first model.nc"])
            destinationText = extractBetween(string(fileread(requestA)), ...
                '"destinations": {','    }');
            testCase.verifyLessThan(strfind(destinationText,'"primary"'),strfind(destinationText,'"secondary"'))
        end

        function writesBarotropicQGGeometryAndIntegrationForms(testCase)
            cases = {
                [5 4], 0, false; ...
                [6 5], 1, true; ...
                [5 6], 1, false; ...
                [6 4], 0, true};
            forms = {
                "fixed-rk4", "explicit"; ...
                "fixed-rk4", "cfl"; ...
                "adaptive-rk23", "adaptive"; ...
                "adaptive-rk45", "adaptive"; ...
                "adaptive-rk78", "adaptive"};
            for iCase = 1:size(cases,1)
                Nxy = cases{iCase,1};
                j = cases{iCase,2};
                shouldAntialias = cases{iCase,3};
                caseName = sprintf("qg-%dx%d-j%d-a%d",Nxy(1),Nxy(2),j,shouldAntialias);
                modelPath = fullfile(testCase.TemporaryFolder,caseName+".nc");
                testCase.createBarotropicQGFixture(modelPath,Nxy,j,shouldAntialias);
                ncwriteatt(modelPath,"/","portableFileIdentifier",caseName);
                for iForm = 1:size(forms,1)
                    method = forms{iForm,1};
                    policy = forms{iForm,2};
                    requestPath = fullfile(testCase.TemporaryFolder,caseName+"-"+method+"-"+policy+".json");
                    if policy == "explicit"
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,initialStep=.0025);
                    elseif policy == "cfl"
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,cfl=.25,timeStepConstraint="advective");
                    else
                        WVModel.writePortableRunRequest(requestPath,modelPath, ...
                            schemaVersion=2,method=method,finalTime=.01,initialStep=.0025,maximumStep=.005, ...
                            relativeTolerance=1e-6,absoluteToleranceScale=1e-9);
                    end
                    document = jsondecode(fileread(requestPath));
                    testCase.verifyEqual(string(document.integration.method),method)
                    testCase.verifyEqual(string(document.modelFiles),modelPath)
                    requestText = lower(string(fileread(requestPath)));
                    for forbidden = ["barotropic","transform","forcing","observer","restart","state"]
                        testCase.verifyFalse(contains(requestText,'"'+forbidden+'"'))
                    end
                end
            end
        end

        function barotropicQGSiblingsGroupsPoliciesAndMappings(testCase)
            primary = fullfile(testCase.TemporaryFolder,"qg primary.nc");
            secondary = fullfile(testCase.TemporaryFolder,"qg secondary.nc");
            testCase.createBarotropicQGFixture(primary,[7 6],0,false,shouldAddSecondGroup=true);
            copyfile(primary,secondary)
            ncwriteatt(primary,"/","portableFileIdentifier","primary");
            ncwriteatt(secondary,"/","portableFileIdentifier","secondary");

            createMap = configureDictionary("string","string");
            createMap(["secondary","primary"]) = ["created-secondary.nc","created-primary.nc"];
            createRequest = fullfile(testCase.TemporaryFolder,"qg-create.json");
            WVModel.writePortableRunRequest(createRequest,[secondary primary], ...
                schemaVersion=2,method="fixed-rk4",finalTime=.01,initialStep=.0025, ...
                outputPolicy="create",destinations=createMap);
            document = jsondecode(fileread(createRequest));
            testCase.verifyEqual(string(document.modelFiles),[secondary;primary])

            for policy = ["replace","append"]
                firstDestination = fullfile(testCase.TemporaryFolder,policy+"-primary.nc");
                secondDestination = fullfile(testCase.TemporaryFolder,policy+"-secondary.nc");
                copyfile(primary,firstDestination)
                copyfile(secondary,secondDestination)
                destinationMap = configureDictionary("string","string");
                destinationMap(["primary","secondary"]) = [firstDestination,secondDestination];
                requestPath = fullfile(testCase.TemporaryFolder,"qg-"+policy+".json");
                WVModel.writePortableRunRequest(requestPath,[primary secondary], ...
                    schemaVersion=2,method="adaptive-rk45",finalTime=.01, ...
                    initialStep=.0025,maximumStep=.005,relativeTolerance=1e-6,absoluteToleranceScale=1e-9, ...
                    outputPolicy=policy,destinations=destinationMap);
                testCase.verifyEqual(string(jsondecode(fileread(requestPath)).output.policy),policy)
            end
        end

        function rejectsInvalidBarotropicQGMetadataAndStateOwnership(testCase)
            request = fullfile(testCase.TemporaryFolder,"qg-invalid.json");
            missing = fullfile(testCase.TemporaryFolder,"qg-missing.nc");
            fieldsOnly = fullfile(testCase.TemporaryFolder,"qg-fields-only.nc");
            testCase.createBarotropicQGFixture(missing,[6 5],1,false,shouldUseLinearDynamics=true);
            testCase.createBarotropicQGFixture(fieldsOnly,[6 5],1,false, ...
                shouldUseLinearDynamics=true,shouldAddA0Field=true);
            for path = [missing fieldsOnly]
                testCase.verifyError(@()testCase.writeFixedQGRequest(request,path), ...
                    "WaveVortexModel:PortableRunRequestContract")
            end

            ambiguous = fullfile(testCase.TemporaryFolder,"qg-ambiguous.nc");
            copyfile(testCase.BarotropicQGPath,ambiguous)
            compactLength = ncinfo(ambiguous,"kl").Size;
            nccreate(ambiguous,"/wave-vortex/A0",Dimensions={"kl",compactLength},Datatype="double");
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,ambiguous), ...
                "WaveVortexModel:PortableRunRequestContract")

            dummy = fullfile(testCase.TemporaryFolder,"qg-dummy-ap.nc");
            copyfile(testCase.BarotropicQGPath,dummy)
            nccreate(dummy,"/wave-vortex/Ap",Dimensions={"kl",compactLength},Datatype="double");
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,dummy), ...
                "WaveVortexModel:PortableRunRequestContract")

            invalidGeometry = fullfile(testCase.TemporaryFolder,"qg-invalid-geometry.nc");
            copyfile(testCase.BarotropicQGPath,invalidGeometry)
            x = ncread(invalidGeometry,"x");
            x(2) = x(1);
            ncwrite(invalidGeometry,"x",x)
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,invalidGeometry), ...
                "WaveVortexModel:PortableRunRequestContract")

            invalidJ = fullfile(testCase.TemporaryFolder,"qg-invalid-j.nc");
            copyfile(testCase.BarotropicQGPath,invalidJ)
            ncwrite(invalidJ,"j",2)
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,invalidJ), ...
                "WaveVortexModel:PortableRunRequestContract")

            invalidAntialias = fullfile(testCase.TemporaryFolder,"qg-invalid-antialias.nc");
            copyfile(testCase.BarotropicQGPath,invalidAntialias)
            ncwrite(invalidAntialias,"shouldAntialias",uint8(2))
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,invalidAntialias), ...
                "WaveVortexModel:PortableRunRequestContract")

            invalidVersion = fullfile(testCase.TemporaryFolder,"qg-invalid-version.nc");
            copyfile(testCase.BarotropicQGPath,invalidVersion)
            ncwriteatt(invalidVersion,"/","model_version","5.0.0")
            testCase.verifyError(@()testCase.writeFixedQGRequest(request,invalidVersion), ...
                "WaveVortexModel:PortableRunRequestContract")

            inconsistent = fullfile(testCase.TemporaryFolder,"qg-inconsistent.nc");
            copyfile(testCase.BarotropicQGPath,inconsistent)
            ncwriteatt(inconsistent,"/","portableFileIdentifier","qg-secondary");
            ncwrite(inconsistent,"h",2*ncread(inconsistent,"h"))
            testCase.verifyError(@()WVModel.writePortableRunRequest(request, ...
                [testCase.BarotropicQGPath inconsistent],schemaVersion=2,method="fixed-rk4", ...
                finalTime=.01,initialStep=.0025),"WaveVortexModel:PortableRunRequestInconsistentBundle")
        end

        function rejectsInvalidContractsOptionsAndAliases(testCase)
            request = fullfile(testCase.TemporaryFolder,"request.json");
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=1,method="adaptive-rk45",finalTime=1,initialStep=.1,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6), ...
                "WaveVortexModel:PortableRunRequestIntegration")
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1,cfl=.2), ...
                "WaveVortexModel:PortableRunRequestIntegration")
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="adaptive-rk45",finalTime=1,initialStep=.1,maximumStep=1,relativeTolerance=1e-3,absoluteToleranceScale=1e-6,cfl=.2), ...
                "WaveVortexModel:PortableRunRequestIntegration")
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,[testCase.ModelPath testCase.ModelPath], ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1), ...
                "WaveVortexModel:PortableRunRequestDuplicateFile")
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,[testCase.ModelPath fullfile(testCase.TemporaryFolder,"missing.nc")], ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1), ...
                "WaveVortexModel:PortableRunRequestMissingFile")

            second = fullfile(testCase.TemporaryFolder,"inconsistent.nc");
            copyfile(testCase.ModelPath,second);
            ncwriteatt(second,"/","portableFileIdentifier","secondary");
            ncwrite(second,"Lx",2*ncread(second,"Lx"));
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,[testCase.ModelPath second], ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1), ...
                "WaveVortexModel:PortableRunRequestInconsistentBundle")

            unsupported = fullfile(testCase.TemporaryFolder,"unsupported.nc");
            copyfile(testCase.ModelPath,unsupported);
            ncwriteatt(unsupported,"/","model_version","5.0.0");
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,unsupported, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1), ...
                "WaveVortexModel:PortableRunRequestContract")

            emptyMap = configureDictionary("string","string");
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1,outputPolicy="create",destinations=emptyMap), ...
                "WaveVortexModel:PortableRunRequestDestinations")
            aliasMap = configureDictionary("string","string");
            aliasMap("primary") = testCase.ModelPath;
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1,outputPolicy="replace",destinations=aliasMap), ...
                "WaveVortexModel:PortableRunRequestAlias")
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1,reportPath=request), ...
                "WaveVortexModel:PortableRunRequestAlias")

            configuration = testCase.baseConfiguration();
            configuration.unknownField = true;
            testCase.verifyError(@()WVModel.writePortableRunRequestForTesting( ...
                request,testCase.ModelPath,configuration,""), ...
                "WaveVortexModel:PortableRunRequestOptions")

            appendTarget = fullfile(testCase.TemporaryFolder,"append target.nc");
            copyfile(testCase.ModelPath,appendTarget);
            ncwrite(appendTarget,"Lx",2*ncread(appendTarget,"Lx"));
            appendMap = configureDictionary("string","string");
            appendMap("primary") = appendTarget;
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,testCase.ModelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1, ...
                outputPolicy="append",destinations=appendMap), ...
                "WaveVortexModel:PortableRunRequestInconsistentBundle")

            multiSecond = fullfile(testCase.TemporaryFolder,"multi second.nc");
            copyfile(testCase.ModelPath,multiSecond);
            ncwriteatt(multiSecond,"/","portableFileIdentifier","secondary");
            collisionMap = configureDictionary("string","string");
            collisionMap(["primary","secondary"]) = ["same.nc","./same.nc"];
            testCase.verifyError(@()WVModel.writePortableRunRequest(request,[testCase.ModelPath multiSecond], ...
                schemaVersion=2,method="fixed-rk4",finalTime=1,initialStep=.1, ...
                outputPolicy="create",destinations=collisionMap), ...
                "WaveVortexModel:PortableRunRequestAlias")
        end

        function everyInjectedFailurePreservesOldRequest(testCase)
            request = fullfile(testCase.TemporaryFolder,"request.json");
            oldBytes = uint8('old request bytes')';
            testCase.writeBytes(request,oldBytes);
            configuration = testCase.baseConfiguration();
            beforeFiles = string({dir(testCase.TemporaryFolder).name});
            for stage = ["validation","encoding","temporary-file","replacement","cleanup"]
                testCase.verifyError(@()WVModel.writePortableRunRequestForTesting( ...
                    request,testCase.ModelPath,configuration,stage), ...
                    "WaveVortexModel:PortableRunRequestInjectedFailure")
                testCase.verifyEqual(testCase.fileBytes(request),oldBytes)
                testCase.verifyEqual(string({dir(testCase.TemporaryFolder).name}),beforeFiles)
            end
        end
    end

    methods (Access=private)
        function request = writeRequest(testCase,folder,name,schemaVersion,method,options)
            arguments
                testCase
                folder
                name
                schemaVersion
                method
                options.initialStep = NaN
                options.cfl = NaN
                options.timeStepConstraint = ""
                options.maximumStep = NaN
                options.relativeTolerance = NaN
                options.absoluteToleranceScale = NaN
            end
            request = fullfile(folder,name+".json");
            WVModel.writePortableRunRequest(request,"model bundle Ω.nc", ...
                schemaVersion=schemaVersion,method=method,finalTime=12, ...
                initialStep=options.initialStep,cfl=options.cfl,timeStepConstraint=options.timeStepConstraint, ...
                maximumStep=options.maximumStep,relativeTolerance=options.relativeTolerance, ...
                absoluteToleranceScale=options.absoluteToleranceScale,reportPath=name+" report Ω.json");
            firstBytes = testCase.fileBytes(request);
            WVModel.writePortableRunRequest(request,"model bundle Ω.nc", ...
                schemaVersion=schemaVersion,method=method,finalTime=12, ...
                initialStep=options.initialStep,cfl=options.cfl,timeStepConstraint=options.timeStepConstraint, ...
                maximumStep=options.maximumStep,relativeTolerance=options.relativeTolerance, ...
                absoluteToleranceScale=options.absoluteToleranceScale,reportPath=name+" report Ω.json");
            testCase.verifyEqual(testCase.fileBytes(request),firstBytes)
        end

        function configuration = baseConfiguration(~)
            configuration = struct( ...
                "schemaVersion",2,"method","fixed-rk4","finalTime",12, ...
                "initialStep",.25,"cfl",NaN,"timeStepConstraint","", ...
                "maximumStep",NaN,"relativeTolerance",NaN,"absoluteToleranceScale",NaN, ...
                "outputPolicy","append","destinations",configureDictionary("string","string"), ...
                "fftProvider","reference","threads",1,"reportPath","report.json");
        end

        function path = fixturePath(testCase)
            path = fullfile(testCase.RepositoryRoot,"PortableRuntime","tests","fixtures","time-series-hydrostatic.nc");
        end

        function bytes = fileBytes(~,path)
            fileIdentifier = fopen(path,"rb");
            cleanup = onCleanup(@()fclose(fileIdentifier));
            bytes = fread(fileIdentifier,Inf,"*uint8");
            clear cleanup
        end

        function writeBytes(~,path,bytes)
            fileIdentifier = fopen(path,"wb");
            cleanup = onCleanup(@()fclose(fileIdentifier));
            fwrite(fileIdentifier,bytes,"uint8");
            clear cleanup
        end

        function writeFixedQGRequest(~,requestPath,modelPath)
            WVModel.writePortableRunRequest(requestPath,modelPath, ...
                schemaVersion=2,method="fixed-rk4",finalTime=.01,initialStep=.0025);
        end

        function createBarotropicQGFixture(~,path,Nxy,j,shouldAntialias,options)
            arguments
                ~
                path (1,1) string
                Nxy (1,2) double
                j (1,1) double
                shouldAntialias (1,1) logical
                options.shouldUseLinearDynamics (1,1) logical = false
                options.shouldAddA0Field (1,1) logical = false
                options.shouldAddSecondGroup (1,1) logical = false
            end
            wvt = WVTransformBarotropicQG([15000 9000],Nxy,h=.8,j=j, ...
                g=9.80665,planetaryRadius=6.3712e6,rotationRate=7.292115e-5, ...
                latitude=33,shouldAntialias=shouldAntialias);
            model = WVModel(wvt,shouldUseLinearDynamics=options.shouldUseLinearDynamics);
            if options.shouldAddA0Field
                model.eulerianObservingSystem.addNetCDFOutputVariables('A0');
            end
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=.005,shouldOverwriteExisting=true);
            if options.shouldAddSecondGroup
                group = outputFile.addNewEvenlySpacedOutputGroup("fields",outputInterval=.005);
                group.addObservingSystem(WVEulerianFields(model,fieldNames={'u','qgpv'}));
            end
            model.setupIntegrator(integratorType="fixed",deltaT=.0025);
            model.integrateToTime(.005,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
        end
    end
end
