classdef TestBenchmarkWebsiteDocumentation < matlab.unittest.TestCase
    properties
        repositoryRoot
        temporaryFolder
    end

    methods (TestClassSetup)
        function addToolsPath(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(testCase.repositoryRoot,"tools")));
        end
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function emptyCatalogBuildsUsefulPageWithoutAssets(testCase)
            [root,buildFolder] = testCase.createFixture("empty");
            testCase.writeCatalog(root,struct([]));

            generateBenchmarkWebsiteDocumentation(root,buildFolder);

            page = string(fileread(fullfile(buildFolder,"benchmarks.md")));
            testCase.verifySubstring(page,"No approved scaling datasets have been published yet");
            testCase.verifySubstring(page,"No approved computer results have been published yet");
            testCase.verifySubstring(page,"No approved result files have been published yet");
            testCase.verifySubstring(page,"No approved matched three-interface result has been published yet");
            testCase.verifyFalse(isfolder(fullfile(buildFolder,"assets","benchmarks")));
            testCase.verifyFalse(contains(page,"Performance across releases"));
        end

        function matlabCppAndMissingCasesRenderDeterministically(testCase)
            [root,firstBuild] = testCase.createFixture("multi-platform");
            cases = standardCases(1);
            first = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260810T120000Z","matlab","builtin","m5-max","M5 Max","4.2.1","2026-08-10T12:00:00Z",cases);
            second = publishedDataset("scaling-standard-v1--matlab-builtin--zen4--20260810T130000Z","matlab","builtin","zen4","Zen 4 workstation","4.2.1","2026-08-10T13:00:00Z",standardCases(1.4));
            cppCases = standardCases(0.6);
            cppCases{1} = unavailableCase(cppCases{1},"Transform is not implemented by this backend.");
            cppCases(2) = [];
            third = publishedDataset("scaling-standard-v1--cpp-fftw--m5-max--20260810T140000Z","cpp","fftw","m5-max","M5 Max","0.1.0","2026-08-10T14:00:00Z",cppCases);
            entries = testCase.publishDatasets(root,{first,second,third});
            testCase.writeCatalog(root,entries);

            generateBenchmarkWebsiteDocumentation(root,firstBuild);
            secondBuild = fullfile(root,"second-build");
            mkdir(secondBuild);
            copyfile(fullfile(testCase.repositoryRoot,"Documentation","WebsiteDocumentation","benchmarks.md"),fullfile(secondBuild,"benchmarks.md"));
            generateBenchmarkWebsiteDocumentation(root,secondBuild);

            comparison = compareDocumentationTrees(firstBuild,secondBuild);
            testCase.verifyTrue(comparison.IsEqual,strjoin(comparison.Substantive,newline));
            page = string(fileread(fullfile(firstBuild,"benchmarks.md")));
            testCase.verifySubstring(page,"M5 Max");
            testCase.verifySubstring(page,"Zen 4 workstation");
            testCase.verifySubstring(page,"WaveVortexModel C++");
            testCase.verifySubstring(page,"Transform is not implemented");
            testCase.verifySubstring(page,"Case was not recorded in this dataset");
            testCase.verifySubstring(page,"<table>");
            testCase.verifySubstring(page,"<th scope=""col"">");
            testCase.verifyFalse(contains(page,"| --- |"));
            testCase.verifyFalse(contains(page,"Performance across releases"));

            chartNames = ["runtime-horizontal" "runtime-vertical" "memory-horizontal" "memory-vertical"];
            for chartName = chartNames
                chartPath = fullfile(firstBuild,"assets","benchmarks",chartName + ".svg");
                testCase.verifyTrue(isfile(chartPath));
                chart = string(fileread(chartPath));
                testCase.verifySubstring(chart,"<title");
                testCase.verifySubstring(chart,"<desc");
                testCase.verifySubstring(chart,"viewBox");
                testCase.verifySubstring(chart,"constant nonhydrostatic transform");
                testCase.verifyEqual(numel(strfind(chart,'class="legend-label"')),3);
                testCase.verifySubstring(chart,'class="legend-label" x="72" y="485"');
            end
            horizontalChart = string(fileread(fullfile(firstBuild,"assets","benchmarks","runtime-horizontal.svg")));
            testCase.verifySubstring(horizontalChart,"Horizontal grid size (Nx = Ny)");
            testCase.verifySubstring(horizontalChart,"Median nonlinear-flux evaluation time (s)");
            verticalChart = string(fileread(fullfile(firstBuild,"assets","benchmarks","runtime-vertical.svg")));
            testCase.verifySubstring(verticalChart,"Vertical grid size (Nz)");
            for dataset = [first second third]
                testCase.verifyTrue(isfile(fullfile(firstBuild,"benchmarks","data",dataset.datasetId + ".json")));
                testCase.verifyTrue(isfile(fullfile(firstBuild,"benchmarks","raw",dataset.datasetId + ".json")));
            end
            report = validateWebsiteDocumentation(firstBuild,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            testCase.verifyTrue(report.IsValid,strjoin(report.Diagnostics,newline));
        end

        function comparableVersionsProduceReleaseHistory(testCase)
            [root,buildFolder] = testCase.createFixture("history");
            first = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260810T120000Z","matlab","builtin","m5-max","M5 Max","4.2.1","2026-08-10T12:00:00Z",standardCases(1));
            second = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260811T120000Z","matlab","builtin","m5-max","M5 Max","4.2.2","2026-08-11T12:00:00Z",standardCases(0.9));
            entries = testCase.publishDatasets(root,{first,second});
            testCase.writeCatalog(root,entries);

            generateBenchmarkWebsiteDocumentation(root,buildFolder);

            page = string(fileread(fullfile(buildFolder,"benchmarks.md")));
            testCase.verifySubstring(page,"## Performance across releases");
            testCase.verifySubstring(page,"4.2.1");
            testCase.verifySubstring(page,"4.2.2");
            testCase.verifyFalse(contains(page,"BENCHMARKS:AT_GLANCE"));
        end

        function matchedInterfaceResultUsesConsistentPublicPresentation(testCase)
            [root,buildFolder] = testCase.createFixture("interfaces");
            [firstEntry,first] = testCase.publishInterfaceDataset(root,[256 256 129],"20260815T120000Z");
            [secondEntry,second] = testCase.publishInterfaceDataset(root,[512 512 257],"20260815T130000Z");
            second.provider.moduleSHA256 = repmat('e',1,64);
            first = compositeInterfaceFixture(first);
            second = compositeInterfaceFixture(second);
            testCase.writeJson(fullfile(root,firstEntry.artifact),first);
            testCase.writeJson(fullfile(root,secondEntry.artifact),second);
            testCase.writeCatalog(root,struct([]),[firstEntry secondEntry]);

            generateBenchmarkWebsiteDocumentation(root,buildFolder);

            page = string(fileread(fullfile(buildFolder,"benchmarks.md")));
            comparison = extractBetween(page,"<!-- BENCHMARKS:INTERFACE_COMPARISON:START -->","<!-- BENCHMARKS:INTERFACE_COMPARISON:END -->");
            testCase.verifySubstring(comparison,"256×256×129");
            testCase.verifySubstring(comparison,"512×512×257");
            testCase.verifySubstring(page,"Adaptive RK3(2) + output");
            testCase.verifySubstring(comparison,"2.000×");
            testCase.verifySubstring(comparison,"−25.0% memory");
            testCase.verifySubstring(comparison,"MATLAB + compiled core");
            testCase.verifySubstring(comparison,"Standalone C++");
            testCase.verifyFalse(contains(comparison,"Process wall"));
            testCase.verifyFalse(contains(comparison,"Maximum error"));
            testCase.verifyFalse(contains(comparison,"incremental"));
            testCase.verifyFalse(contains(comparison,"final RSS","IgnoreCase",true));
            testCase.verifyFalse(contains(comparison,"PREVIEW-AVAILABLE"));
            testCase.verifySubstring(page,string(first.schemaVersion));
            testCase.verifySubstring(page,"/benchmarks/data/"+first.datasetId+".json");
            testCase.verifySubstring(page,"/benchmarks/data/"+second.datasetId+".json");
            testCase.verifySubstring(page,"External archives:");
            testCase.verifyTrue(isfile(fullfile(buildFolder,"benchmarks","data",first.datasetId+".json")));
            testCase.verifyTrue(isfile(fullfile(buildFolder,"benchmarks","data",second.datasetId+".json")));
            testCase.verifyFalse(isfile(fullfile(buildFolder,"benchmarks","raw",first.datasetId+".json")));
        end

        function integratorStudyRendersOneRuntimeAndRSSMatrix(testCase)
            [root,buildFolder] = testCase.createFixture("integrator-study");
            [firstEntry,first] = testCase.publishIntegratorStudyDataset(root,[256 256 129],"20260825T120000Z");
            testCase.writeCatalog(root,struct([]),firstEntry);

            generateBenchmarkWebsiteDocumentation(root,buildFolder);

            page = string(fileread(fullfile(buildFolder,"benchmarks.md")));
            comparison = extractBetween(page,"<!-- BENCHMARKS:INTERFACE_COMPARISON:START -->","<!-- BENCHMARKS:INTERFACE_COMPARISON:END -->");
            testCase.verifyEqual(numel(strfind(comparison,"256×256×129")),8)
            testCase.verifyFalse(contains(comparison,"512×512×257"))
            testCase.verifyFalse(contains(comparison,">Hydrostatic<"))
            testCase.verifySubstring(comparison,"Nonhydrostatic")
            testCase.verifySubstring(comparison,"Fixed RK4")
            testCase.verifySubstring(comparison,"ode23 / RK3(2)")
            testCase.verifySubstring(comparison,"ode45 / RK5(4)")
            testCase.verifySubstring(comparison,"ode78 / RK8(7)")
            testCase.verifySubstring(comparison,"Coefficients · endpoint only")
            testCase.verifySubstring(comparison,"Composite graph · interior dense output")
            testCase.verifySubstring(comparison,"1 s / 1 GiB")
            testCase.verifyFalse(contains(comparison,"×)"))
            testCase.verifyFalse(contains(comparison,"incremental","IgnoreCase",true))
            testCase.verifyFalse(contains(comparison,"process wall","IgnoreCase",true))
            testCase.verifySubstring(comparison,"exact standalone workspace ledgers")
            testCase.verifyTrue(isfile(fullfile(buildFolder,"benchmarks","data",first.datasetId+".json")))
        end

        function incompleteOrIncompatibleInterfacePairFails(testCase)
            [root,buildFolder] = testCase.createFixture("incomplete-interfaces");
            [firstEntry,~] = testCase.publishInterfaceDataset(root,[256 256 129],"20260815T120000Z");
            testCase.writeCatalog(root,struct([]),firstEntry);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:IncompleteInterfaceComparison");

            [secondEntry,second] = testCase.publishInterfaceDataset(root,[512 512 257],"20260815T130000Z");
            second.cases{2}.contract.deltaT = 2e-3;
            testCase.writeJson(fullfile(root,secondEntry.artifact),second);
            testCase.writeCatalog(root,struct([]),[firstEntry secondEntry]);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:IncompleteInterfaceComparison");
        end

        function missingSuiteDatasetIsExplicitlyUnavailable(testCase)
            [root,buildFolder] = testCase.createFixture("missing-suite");
            standard = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260810T120000Z","matlab","builtin","m5-max","M5 Max","4.2.1","2026-08-10T12:00:00Z",standardCases(1));
            large = publishedDataset("scaling-large-v1--matlab-builtin--zen4--20260810T130000Z","matlab","builtin","zen4","Zen 4 workstation","4.2.1","2026-08-10T13:00:00Z",standardCases(1.4));
            large.benchmark.suiteId = "scaling-large-v1";
            entries = testCase.publishDatasets(root,{standard,large});
            testCase.writeCatalog(root,entries);

            generateBenchmarkWebsiteDocumentation(root,buildFolder);

            page = string(fileread(fullfile(buildFolder,"benchmarks.md")));
            testCase.verifySubstring(page,"No scaling-large-v1 dataset was collected for this environment.");
            testCase.verifySubstring(page,"No scaling-standard-v1 dataset was collected for this environment.");
            testCase.verifySubstring(page,"MATLAB R2026a Update 4");
            testCase.verifySubstring(page,"<td>scaling-large-v1</td>");
            testCase.verifySubstring(page,"<td>Constant nonhydrostatic</td>");
            testCase.verifySubstring(page,"<td>scaling-standard-v1</td>");
            chart = string(fileread(fullfile(buildFolder,"assets","benchmarks","runtime-horizontal.svg")));
            testCase.verifySubstring(chart,"M5 Max — MATLAB R2026a Update 4 — standard");
            testCase.verifySubstring(chart,"Zen 4 workstation — MATLAB R2026a Update 4 — large");
            testCase.verifyFalse(contains(chart,"Barotropic QG"));
        end

        function unsafeMissingDuplicateAndMismatchedArtifactsFail(testCase)
            [root,buildFolder] = testCase.createFixture("invalid");
            dataset = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260810T120000Z","matlab","builtin","m5-max","M5 Max","4.2.1","2026-08-10T12:00:00Z",standardCases(1));
            entries = testCase.publishDatasets(root,{dataset});

            unsafe = entries;
            unsafe.artifact = "../outside.json";
            testCase.writeCatalog(root,unsafe);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:UnsafeBenchmarkPath");

            missing = entries;
            missing.artifact = "Benchmarks/results/published/missing.json";
            testCase.writeCatalog(root,missing);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:MissingBenchmarkArtifact");

            mismatch = entries;
            mismatch.datasetId = "scaling-standard-v1--matlab-builtin--other--20260810T120000Z";
            testCase.writeCatalog(root,mismatch);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:PublishedBenchmarkIdentityMismatch");

            duplicate = [entries entries];
            testCase.writeCatalog(root,duplicate);
            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:DuplicatePublishedBenchmark");
        end

        function conflictingComparableCasesFail(testCase)
            [root,buildFolder] = testCase.createFixture("conflict");
            first = publishedDataset("scaling-standard-v1--matlab-builtin--m5-max--20260810T120000Z","matlab","builtin","m5-max","M5 Max","4.2.1","2026-08-10T12:00:00Z",standardCases(1));
            conflictingCases = standardCases(1.2);
            conflictingCases{1}.configuration.seed = 999;
            second = publishedDataset("scaling-standard-v1--matlab-builtin--zen4--20260810T130000Z","matlab","builtin","zen4","Zen 4 workstation","4.2.1","2026-08-10T13:00:00Z",conflictingCases);
            entries = testCase.publishDatasets(root,{first,second});
            testCase.writeCatalog(root,entries);

            testCase.verifyError(@()generateBenchmarkWebsiteDocumentation(root,buildFolder),"WaveVortexModel:ConflictingPublishedBenchmarkCase");
        end
    end

    methods (Access=private)
        function [root,buildFolder] = createFixture(testCase,name)
            root = fullfile(testCase.temporaryFolder,name);
            buildFolder = fullfile(root,"build");
            mkdir(fullfile(root,"Benchmarks","results"));
            mkdir(buildFolder);
            copyfile(fullfile(testCase.repositoryRoot,"Documentation","WebsiteDocumentation","benchmarks.md"),fullfile(buildFolder,"benchmarks.md"));
        end

        function entries = publishDatasets(testCase,root,datasets)
            entries = repmat(struct("datasetId","","artifact",""),1,numel(datasets));
            publishedFolder = fullfile(root,"Benchmarks","results","published");
            rawFolder = fullfile(root,"Benchmarks","results","raw");
            mkdir(publishedFolder);
            mkdir(rawFolder);
            for iDataset = 1:numel(datasets)
                dataset = datasets{iDataset};
                artifact = "Benchmarks/results/published/" + dataset.datasetId + ".json";
                rawArtifact = "Benchmarks/results/raw/" + dataset.datasetId + ".json";
                dataset.provenance.rawArtifact = rawArtifact;
                testCase.writeJson(fullfile(root,artifact),dataset);
                testCase.writeJson(fullfile(root,rawArtifact),struct("fixture",true));
                entries(iDataset) = struct("datasetId",dataset.datasetId,"artifact",artifact);
            end
        end

        function [entry,dataset] = publishInterfaceDataset(testCase,root,resolution,timestamp)
            interfaces = { ...
                struct("id","matlab-builtin","processWallSeconds",4,"interfaceTotalSeconds",3,"integrationSeconds",2,"totalPeakRSSBytes",4*2^30,"incrementalPeakRSSBytes",2*2^30,"processWallRatio",1,"integrationRatio",1,"totalRSSRatio",1,"incrementalRSSRatio",1,"processWallSamplesSeconds",[3.9 4 4.1],"integrationSamplesSeconds",[1.9 2 2.1]), ...
                struct("id","matlab-compiled","processWallSeconds",3,"interfaceTotalSeconds",2,"integrationSeconds",1,"totalPeakRSSBytes",3*2^30,"incrementalPeakRSSBytes",1*2^30,"processWallRatio",0.75,"integrationRatio",0.5,"totalRSSRatio",0.75,"incrementalRSSRatio",0.5,"processWallSamplesSeconds",[2.9 3 3.1],"integrationSamplesSeconds",[0.9 1 1.1]), ...
                struct("id","standalone-compiled","processWallSeconds",1,"interfaceTotalSeconds",0.8,"integrationSeconds",0.5,"totalPeakRSSBytes",1*2^30,"incrementalPeakRSSBytes",0.5*2^30,"processWallRatio",0.25,"integrationRatio",0.25,"totalRSSRatio",0.25,"incrementalRSSRatio",0.25,"processWallSamplesSeconds",[0.9 1 1.1],"integrationSamplesSeconds",[0.4 0.5 0.6])};
            contract = struct("Nxyz",resolution,"Lxyz",[15000 15000 1300],"forcing","default WVNonlinearAdvection","shouldAntialias",true,"integrator","none","deltaT",1e-3,"finalTime",2e-3,"relativeTolerance",1e-3,"absoluteTolerance",1e-6,"outputInterval",5e-4,"observerGraph","fixture","processRunCount",3,"warmupCount",0,"samplesPerProcess",1);
            definitions = ["nonlinear-flux" "fixed-rk4-continuation" "adaptive-rk23-observer-output"];
            cases = cell(1,3);
            categories = arrayfun(@(name)struct("name",name,"variableCount",1,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"passed",true),["coefficients" "eulerianFields" "moorings" "particles" "tracers" "times"]);
            graph = struct("passed",true,"variableCount",6,"recordCount",12,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"categories",categories);
            for iCase = 1:3
                cases{iCase} = struct("id",definitions(iCase),"operation","model-continuation","contract",contract,"interfaces",{interfaces},"correctness",struct("passed",true,"maximumRelativeError",1e-14,"outputAgreementPassed",true,"completeOutputGraph",graph));
            end
            datasetId = "three-interface--lyra--"+timestamp;
            provenance = struct("rawSchemaVersion","three-interface-benchmark-v1","externalArchive",struct("fileName",datasetId+".json.gz","sha256",repmat('d',1,64),"compressedBytes",1024));
            source = struct("repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",repmat('a',1,40),"tree",repmat('b',1,40),"sourceDirty",false,"version","unreleased-preview");
            platform = struct("id","lyra","displayName","Lyra (Apple M4 Max)","processor","Apple M4 Max","physicalMemoryBytes",128*2^30,"os","macOS","architecture","maca64","matlabVersion","R2025b Update 4","threadCount",16);
            provider = struct("id","native-neon-pthreads","version","3.3.11","threadBackend","pthreads","moduleSHA256",repmat('c',1,64),"identityValidated",true,"openMPDetected",false);
            dataset = struct("schemaVersion","published-three-interface-v1","datasetId",datasetId,"collectedAt","2026-08-15T12:00:00Z","source",source,"platform",platform,"provider",provider,"provenance",provenance,"cases",{cases});
            artifact = "Benchmarks/results/published/"+datasetId+".json";
            testCase.writeJson(fullfile(root,artifact),dataset);
            entry = struct("datasetId",datasetId,"artifact",artifact);
        end

        function [entry,dataset] = publishIntegratorStudyDataset(testCase,root,resolution,timestamp)
            interfaceIds = ["matlab-builtin" "matlab-compiled" "standalone-compiled"];
            interfaces = cell(1,3);
            for iInterface = 1:3
                diagnostics = struct("controls",struct("requested","fixed-rk4","actual","fixed-rk4","matched",true),"methodWork",struct("acceptedStepCounts",10,"rejectedStepCounts",0,"rhsEvaluationCounts",40,"denseOutputEvaluationCounts",0),"integratorStorage",struct("exact",iInterface==3),"stateSizedBuffers",struct([]),"memory",struct("boundary","integration-phase-total-live-process-tree-rss"));
                interfaces{iInterface} = struct("id",interfaceIds(iInterface),"providerId",conditional(iInterface==1,"matlab-builtin","native-neon-pthreads"),"integrationSeconds",1/iInterface,"totalPeakRSSBytes",iInterface*2^30,"integrationRatio",1/iInterface,"totalRSSRatio",iInterface,"integrationSamplesSeconds",[0.9 1 1.1]/iInterface,"totalPeakRSSSamplesBytes",iInterface*[0.9 1 1.1]*2^30,"diagnostics",diagnostics);
            end
            categories = arrayfun(@(name)struct("name",name,"variableCount",1,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"passed",true),["coefficients" "eulerianFields" "moorings" "particles" "tracers" "times"]);
            graph = struct("passed",true,"variableCount",6,"recordCount",12,"maximumAbsoluteError",0,"maximumRelativeError",1e-14,"categories",categories);
            cases = cell(1,8);
            iCase = 0;
            physicalConfiguration = "nonhydrostatic";
            for integrator = ["fixed-rk4" "adaptive-rk23" "adaptive-rk45" "adaptive-rk78"]
                for workload = ["coefficient-endpoint" "composite-dense-output"]
                    iCase = iCase+1;
                    identifier = physicalConfiguration+"--"+integrator+"--"+workload;
                    contract = struct("Nxyz",resolution,"Lxyz",[150e3 150e3 1300],"physicalConfiguration",physicalConfiguration,"isHydrostatic",false,"workload",workload,"forcing","default WVNonlinearAdvection","shouldAntialias",true,"integrator",integrator,"deltaT",128,"integrationStepCount",56,"finalTime",7168,"relativeTolerance",1e-3,"absoluteToleranceScale",1e-6,"absoluteToleranceEvidence","component hashes","initialStep",295.7935799274,"maximumStep",716.8,"outputInterval",conditional(workload=="coefficient-endpoint",7168,128),"denseOutputPointsPerStep",4,"observerGraph",workload,"processRunCount",3,"warmupCount",0,"samplesPerProcess",1);
                    correctness = struct("passed",true,"maximumRelativeError",1e-14,"outputAgreementPassed",true,"endpointTrajectoryAgreementPassed",true,"completeOutputGraph",graph);
                    cases{iCase} = struct("id",identifier,"physicalConfiguration",physicalConfiguration,"workload",workload,"integrator",integrator,"contract",contract,"interfaces",{interfaces},"correctness",correctness);
                end
            end
            datasetId = "three-interface--donut--"+timestamp;
            source = struct("repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",repmat('a',1,40),"tree",repmat('b',1,40),"sourceDirty",false,"version","unreleased-preview");
            platform = struct("id","donut","displayName","Donut (Apple M5 Max)","processor","Apple M5 Max","physicalMemoryBytes",48*2^30,"os","macOS","architecture","maca64","matlabVersion","R2026a Update 4","threadCount",18);
            provider = struct("id","native-neon-pthreads","version","3.3.11","threadBackend","pthreads","scope","compiled-interfaces-only","moduleSHA256",repmat('c',1,64),"identityValidated",true,"openMPDetected",false);
            fixtures = [struct("physicalConfiguration","nonhydrostatic","workload","coefficient-endpoint","sha256",repmat('3',1,64)); struct("physicalConfiguration","nonhydrostatic","workload","composite-dense-output","sha256",repmat('4',1,64))];
            archive = struct("fileName",datasetId+".json.gz","sha256",repmat('d',1,64),"compressedBytes",4096,"location","external sibling archive");
            provenance = struct("rawSchemaVersion","three-interface-benchmark-v2","rawArtifactSHA256",repmat('e',1,64),"externalArchive",archive,"fixtures",fixtures);
            dataset = struct("schemaVersion","published-three-interface-v3","datasetId",datasetId,"collectedAt","2026-08-25T12:00:00Z","studyId","integrator-runtime-memory-v1","source",source,"platform",platform,"provider",provider,"provenance",provenance,"cases",{cases});
            artifact = "Benchmarks/results/published/"+datasetId+".json";
            testCase.writeJson(fullfile(root,artifact),dataset);
            entry = struct("datasetId",datasetId,"artifact",artifact);
        end

        function writeCatalog(testCase,root,publishedDatasets,interfaceComparisons)
            if nargin < 4
                interfaceComparisons = struct([]);
            end
            references = [ ...
                struct("suiteId","core-v1","backendId","builtin","rawArtifact","Benchmarks/results/reference/core.json"), ...
                struct("suiteId","scaling-standard-v1","backendId","builtin","rawArtifact","Benchmarks/results/reference/standard.json"), ...
                struct("suiteId","scaling-large-v1","backendId","builtin","rawArtifact","Benchmarks/results/reference/large.json")];
            catalog = struct("schemaVersion","benchmark-catalog-v1","scoringReferences",references,"publishedDatasets",publishedDatasets,"interfaceComparisons",interfaceComparisons);
            testCase.writeJson(fullfile(root,"Benchmarks","results","catalog.json"),catalog);
        end

        function writeJson(~,path,value)
            parent = fileparts(path);
            if ~isfolder(parent)
                mkdir(parent);
            end
            fileId = fopen(path,"w");
            cleanup = onCleanup(@()fclose(fileId));
            fprintf(fileId,"%s\n",jsonencode(value,PrettyPrint=true));
            clear cleanup
        end
    end
end

function dataset = compositeInterfaceFixture(dataset)
provider = dataset.provider;
evidence = struct("datasetId",dataset.datasetId,"collectedAt",dataset.collectedAt,"source",dataset.source,"provider",provider,"externalArchive",dataset.provenance.externalArchive);
for iCase = 1:numel(dataset.cases)
    dataset.cases{iCase}.evidence = evidence;
end
dataset.schemaVersion = "published-three-interface-v2";
dataset.provider.moduleSHA256 = "per-case-evidence";
dataset.provider.moduleIdentityScope = "case-evidence";
dataset.provenance = struct;
dataset.provenance.composition = "frozen-valid-v1-plus-corrected-adaptive";
dataset.provenance.sourceDatasets = [evidence evidence];
end

function dataset = publishedDataset(datasetId,implementationId,backendId,platformId,platformName,version,collectedAt,cases)
benchmark = struct("suiteId","scaling-standard-v1","suiteVersion",1,"operation","nonlinearFlux","correctnessTolerance",1e-12);
displayName = "WaveVortexModel MATLAB";
toolchain = struct("kind","matlab","name","MATLAB","version","R2026a Update 4","details",struct("release","R2026a"));
if implementationId == "cpp"
    displayName = "WaveVortexModel C++";
    toolchain = struct("kind","cpp","name","Clang","version","18.0.0","details",struct("buildType","Release"));
end
implementation = struct("id",implementationId,"displayName",displayName,"version",version,"repository","https://github.com/JeffreyEarly/wave-vortex-model","commit",repmat('a',1,40),"backend",backendId,"sourceDirty",false);
platform = struct("id",platformId,"displayName",platformName,"processor",platformName + " processor","physicalMemoryBytes",64*2^30,"os","Example OS","architecture","example64","threadCount",16);
provenance = struct("rawArtifact","placeholder.json","rawSchemaVersion","1.1.0");
dataset = struct("schemaVersion","published-benchmark-v1","datasetId",datasetId,"collectedAt",collectedAt,"benchmark",benchmark,"implementation",implementation,"platform",platform,"toolchain",toolchain,"provenance",provenance,"cases",{cases});
end

function cases = standardCases(scale)
definitions = { ...
    "constant-nonhydrostatic-64x64x65","constant-nonhydrostatic",[64 64 65]; ...
    "constant-nonhydrostatic-128x128x65","constant-nonhydrostatic",[128 128 65]; ...
    "constant-nonhydrostatic-256x256x65","constant-nonhydrostatic",[256 256 65]; ...
    "constant-nonhydrostatic-128x128x33","constant-nonhydrostatic",[128 128 33]; ...
    "constant-nonhydrostatic-128x128x129","constant-nonhydrostatic",[128 128 129]; ...
    "barotropic-qg-128x128","barotropic-qg",[128 128]; ...
    "barotropic-qg-256x256","barotropic-qg",[256 256]};
cases = cell(1,size(definitions,1));
for iCase = 1:numel(cases)
    cases{iCase} = completeCase(definitions{iCase,1},definitions{iCase,2},definitions{iCase,3},100+iCase,scale*iCase/100);
end
end

function benchmarkCase = completeCase(id,transformId,Nxyz,seed,medianSeconds)
Lxyz = [15e3 15e3 1300];
if numel(Nxyz) == 2
    Lxyz = Lxyz(1:2);
end
configuration = struct("Lxyz",Lxyz,"Nxyz",Nxyz,"isHydrostatic",transformId ~= "constant-nonhydrostatic","shouldAntialias",true,"seed",seed,"warmupCount",2,"sampleCount",3);
timing = struct("medianSeconds",medianSeconds,"samplesSeconds",medianSeconds*[0.9 1 1.1]);
correctness = struct("passed",true,"relativeError",0);
peak = (2+iScale(Nxyz))*2^30;
memory = struct("status","complete","provider","fixture-rss","baselineProcessBytes",2^30,"peakProcessBytes",peak,"peakIncrementBytes",peak-2^30);
benchmarkCase = struct("id",id,"transformId",transformId,"scoreFamily",transformId,"configuration",configuration,"status","complete","timing",timing,"correctness",correctness,"memory",memory);
end

function scale = iScale(Nxyz)
scale = prod(double(Nxyz))/prod([64 64]);
end

function benchmarkCase = unavailableCase(complete,reason)
benchmarkCase = struct("id",complete.id,"transformId",complete.transformId,"scoreFamily",complete.scoreFamily,"configuration",complete.configuration,"status","unavailable","unavailableReason",reason);
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
