classdef TestPortableFieldSamplingMatrix < matlab.unittest.TestCase
    properties (SetAccess=private)
        root (1,1) string
        folder (1,1) string
        probe (1,1) string
        providers (1,:) string
    end
    methods (TestClassSetup)
        function prepareProbe(testCase)
            testCase.root=string(fileparts(fileparts(mfilename("fullpath"))));
            fixture=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture(PreservingOnFailure=true));
            testCase.folder=string(fixture.Folder);
            stableProbe=string(getenv("WV_STABLE_FORCING_DUMP"));
            testCase.probe=string(getenv("WV_FIELD_SAMPLING_DUMP"));
            if testCase.probe=="" && stableProbe~="", testCase.probe=fullfile(fileparts(stableProbe),"WVStratifiedQGFieldDump"); end
            testCase.providers="reference";
            if getenv("WV_STABLE_FORCING_NATIVE")=="1" || getenv("WV_FIELD_SAMPLING_NATIVE")=="1", testCase.providers=["reference","native"]; end
            if testCase.probe==""
                build=fullfile(testCase.folder,"build");
                [status,output]=cleanSystem("cmake -S "+shellQuote(fullfile(testCase.root,"PortableRuntime"))+" -B "+shellQuote(build)+" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output);
                [status,output]=cleanSystem("cmake --build "+shellQuote(build)+" --parallel 2 --target WVStratifiedQGFieldDump");
                testCase.assertEqual(status,0,output);
                testCase.probe=fullfile(build,"WVStratifiedQGFieldDump");
            end
            testCase.assertTrue(isfile(testCase.probe));
        end
    end
    methods (Test,TestTags="full")
        function declaredPortableSamplingMatchesMatlab(testCase)
            catalog=jsondecode(fileread(fullfile(testCase.root,"PortableRuntime","contracts","portable-variable-catalog-v1.json")));
            families=["constant-hydrostatic","constant-nonhydrostatic","barotropic","stratified-qg","hydrostatic","boussinesq"];
            evidence={};
            for family=families
                for antialias=[false true]
                    configuration=family+"-aa"+double(antialias);
                    rows=catalog.contracts(string({catalog.contracts.configuration})==configuration);
                    positions={}; profiles={}; movingPositions={};
                    for row=reshape(rows,1,[])
                        if string(row.runtimeStatus)~="implemented", continue; end
                        modes=string(row.metadata.samplingModes);
                        if ismember("positions",modes), positions{end+1}=string(row.metadata.name); end %#ok<AGROW>
                        if ismember("fixedVerticalProfiles",modes), profiles{end+1}=string(row.metadata.name); end %#ok<AGROW>
                        if ismember("positions",modes), movingPositions{end+1}=string(row.metadata.name); end %#ok<AGROW>
                    end
                    [wvt,source]=testCase.writeSource(family,antialias,configuration);
                    [x,y,z]=testCase.positionCoordinates(wvt);
                    request=struct(x=x,y=y,z=z,fields={positions});
                    if isempty(positions), request.fields={}; end
                    if isempty(profiles)
                        request.profiles=struct(fields={{}},xIndices=[],yIndices=[]);
                    else
                        profileOrd=0:numel(profiles)-1;
                        request.profiles=struct(fields={profiles},xIndices=1+mod(profileOrd,wvt.Nx),yIndices=1+mod(floor(profileOrd/wvt.Nx),wvt.Ny));
                    end
                    requestPath=fullfile(testCase.folder,configuration+"-request.json"); outputPath=fullfile(testCase.folder,configuration+"-output.json");
                    probeRequest=request;
                    writeJSON(requestPath,probeRequest);
                    for provider=testCase.providers
                        [status,output]=cleanSystem(shellQuote(testCase.probe)+" "+shellQuote(source)+" "+shellQuote(requestPath)+" "+shellQuote(outputPath)+" "+provider);
                        testCase.assertEqual(status,0,configuration+" "+provider+": "+output);
                        actual=jsondecode(fileread(outputPath));
                        for method=["linear","spline"]
                            for name=string(positions)
                                expected=wvt.variableAtPositionWithName(request.x,request.y,request.z,char(name),interpolationMethod=char(method));
                                tolerance=testCase.samplingTolerance(wvt,name,expected);
                                cppGrid=reshape(actual.full.(name),size(wvt.(name)));
                                interpolatedCppGrid=testCase.interpolateFullField(wvt,request.x,request.y,request.z,method,cppGrid);
                                interpolationTolerance=2e-11*max(abs(interpolatedCppGrid(:)))+1e-16;
                                if endsWith(name,"_portable_catalog_forcing"), interpolationTolerance=2e-11*max(abs(cppGrid(:)))+realmin; end
                                for path=["fixed","moving","event"]
                                    if ~isfield(actual.(method).(path),name)
                                        testCase.verifyEqual(path,"moving",configuration+" "+provider+" "+method+" missing "+name+" output");
                                        testCase.verifyFalse(ismember(name,string(movingPositions)),configuration+" "+provider+" "+method+" omitted metadata-declared moving "+name);
                                        continue
                                    end
                                    if path=="moving", testCase.verifyTrue(ismember(name,string(movingPositions)),configuration+" "+provider+" "+method+" unexpectedly emitted moving "+name); end
                                    observed=actual.(method).(path).(name);
                                    testCase.verifyLessThanOrEqual(max(abs(observed(:)-expected(:))),tolerance,configuration+" "+provider+" "+method+" "+path+" "+name);
                                    testCase.verifyLessThanOrEqual(max(abs(observed(:)-interpolatedCppGrid(:))),interpolationTolerance,configuration+" "+provider+" "+method+" "+path+" interpolation of full "+name);
                                    evidence{end+1}=struct(configuration=configuration,provider=provider,method=method,path=path,field=name,maximumError=max(abs(observed(:)-expected(:))),tolerance=tolerance,interpolationMaximumError=max(abs(observed(:)-interpolatedCppGrid(:))),interpolationTolerance=interpolationTolerance); %#ok<AGROW>
                                end
                            end
                            for name=string(profiles)
                                expected=testCase.profileColumn(wvt,name,request.profiles.xIndices(find(string(profiles)==name,1)),request.profiles.yIndices(find(string(profiles)==name,1)));
                                observed=actual.(method).profiles.(name);
                                tolerance=testCase.samplingTolerance(wvt,name,expected);
                                testCase.verifyLessThanOrEqual(max(abs(observed(:)-expected(:))),tolerance,configuration+" "+provider+" "+method+" profile "+name);
                                cppGrid=reshape(actual.full.(name),size(wvt.(name)));
                                fieldIndex=find(string(profiles)==name,1);
                                expectedCppColumn=cppGrid(request.profiles.xIndices(fieldIndex),request.profiles.yIndices(fieldIndex),:);
                                testCase.verifyEqual(observed(:),expectedCppColumn(:),configuration+" "+provider+" profile must select the full-field column for "+name);
                                evidence{end+1}=struct(configuration=configuration,provider=provider,method=method,path="fixedVerticalProfiles",field=name,maximumError=max(abs(observed(:)-expected(:))),tolerance=tolerance,interpolationMaximumError=max(abs(observed(:)-expectedCppColumn(:))),interpolationTolerance=0); %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            testCase.writeEvidence("matrix",evidence);
        end
        function densitySamplingPreservesReferenceSelection(testCase)
            names=["eta_true","ape","apv"];
            evidence={};
            for family=["constant-hydrostatic","constant-nonhydrostatic","hydrostatic","boussinesq"]
                for antialias=[false true]
                    configuration=family+"-aa"+double(antialias);
                    [wvt,source]=testCase.writeSource(family,antialias,configuration);
                    [x,y,z]=testCase.positionCoordinates(wvt);
                    wvt.shouldUseTrueNoMotionProfile=true;
                    actualDisplacement=wvt.eta_true;
                    wvt.shouldUseTrueNoMotionProfile=false;
                    testCase.assertGreaterThan(max(abs(actualDisplacement-wvt.eta_true),[],"all"),1e-8,"Density reference fixture must distinguish actual from initial.");
                    for reference=["actual","initial"]
                        wvt.shouldUseTrueNoMotionProfile=reference=="actual";
                        request=struct(x=x,y=y,z=z,fields=names,reference=reference,profiles=struct(fields=names,xIndices=[1 3 5],yIndices=[1 2 3]));
                        requestPath=fullfile(testCase.folder,configuration+"-density-request.json");
                        outputPath=fullfile(testCase.folder,configuration+"-density-output.json");
                        writeJSON(requestPath,request);
                        for provider=testCase.providers
                            [status,output]=cleanSystem(shellQuote(testCase.probe)+" "+shellQuote(source)+" "+shellQuote(requestPath)+" "+shellQuote(outputPath)+" "+provider);
                            testCase.assertEqual(status,0,configuration+" "+reference+" "+provider+": "+output);
                            actual=jsondecode(fileread(outputPath));
                            for method=["linear","spline"]
                                for index=1:numel(names)
                                    name=names(index);
                                    expected=wvt.variableAtPositionWithName(x,y,z,char(name),interpolationMethod=char(method));
                                    tolerance=testCase.samplingTolerance(wvt,name,expected);
                                    cppGrid=reshape(actual.full.(name),size(wvt.(name)));
                                    interpolatedCppGrid=testCase.interpolateFullField(wvt,x,y,z,method,cppGrid);
                                    for path=["fixed","moving","event"]
                                        observed=actual.(method).(path).(name);
                                        testCase.verifyLessThanOrEqual(max(abs(observed(:)-expected(:))),tolerance,configuration+" "+reference+" "+provider+" "+path+" "+name);
                                        interpolationTolerance=2e-11*max(abs(interpolatedCppGrid(:)))+1e-16;
                                        testCase.verifyLessThanOrEqual(max(abs(observed(:)-interpolatedCppGrid(:))),interpolationTolerance);
                                        evidence{end+1}=struct(configuration=configuration,reference=reference,provider=provider,method=method,path=path,field=name,maximumError=max(abs(observed(:)-expected(:))),tolerance=tolerance,interpolationMaximumError=max(abs(observed(:)-interpolatedCppGrid(:))),interpolationTolerance=interpolationTolerance); %#ok<AGROW>
                                    end
                                    expected=testCase.profileColumn(wvt,name,request.profiles.xIndices(index),request.profiles.yIndices(index));
                                    observed=actual.(method).profiles.(name);
                                    testCase.verifyLessThanOrEqual(max(abs(observed(:)-expected(:))),testCase.samplingTolerance(wvt,name,expected));
                                    cppColumn=cppGrid(request.profiles.xIndices(index),request.profiles.yIndices(index),:);
                                    testCase.verifyEqual(observed(:),cppColumn(:));
                                end
                            end
                        end
                    end
                end
            end
            testCase.writeEvidence("density-reference",evidence);
        end
    end
    methods (Access=private)
        function [wvt,path]=writeSource(testCase,family,antialias,configuration)
            grid=[8 6 9]; if antialias, grid=[7 5 8]; end
            if family=="barotropic"
                wvt=WVTransformBarotropicQG([17000 11000],grid(1:2),j=1,shouldAntialias=antialias);
            elseif family=="stratified-qg"
                wvt=WVTransformStratifiedQG([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
            elseif family=="hydrostatic"
                wvt=WVTransformHydrostatic([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
            elseif family=="boussinesq"
                wvt=WVTransformBoussinesq([17000 11000 1000],grid,Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=antialias);
            else
                wvt=WVTransformConstantStratification([17000 11000 1000],grid,N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=antialias);
            end
            n=reshape(1:numel(wvt.A0),size(wvt.A0));
            if ismember(family,["barotropic","stratified-qg"])
                wvt.A0=1e-6*complex(sin(.7*n),cos(.3*n));
                if isprop(wvt,"K"), wvt.A0(wvt.K==0 & wvt.L==0)=0; end
            else
                wave=wvt.J>0 & (wvt.K.^2+wvt.L.^2)>0; inertial=(wvt.K.^2+wvt.L.^2)==0;
                wvt.Ap=.001*(sin(.7*n)+1i*cos(.3*n))./(1+wvt.J).*(wave|inertial);
                wvt.Am=.002*(sin(.3*n)+1i*cos(.7*n))./(1+wvt.J).*wave;
                wvt.Am(inertial)=conj(wvt.Ap(inertial));
                mda=inertial & wvt.J>0; wvt.A0(mda)=.003*sin(n(mda));
            end
            wvt.t=37; wvt.t0=17; wvt.removeAllForcing();
            held=WVFixedAmplitudeForcing(wvt,name="portable_catalog_forcing");
            held.setGeostrophicForcingCoefficients(wvt.A0);
            if ~ismember(family,["barotropic","stratified-qg"]), held.setWaveForcingCoefficients(wvt.Ap,wvt.Am); end
            wvt.setForcing([WVNonlinearAdvection(wvt),held]);
            wvt.addOperation(SpatialForcingOperation(wvt));
            forcingNames=string(wvt.variableNames);
            forcingNames=forcingNames(endsWith(forcingNames,"_portable_catalog_forcing"));
            forcingMagnitude=0;
            for name=reshape(forcingNames,1,[]), forcingMagnitude=max(forcingMagnitude,max(abs(wvt.(name)),[],"all")); end
            testCase.assertGreaterThan(forcingMagnitude,0,"The forcing sampling fixture must have a nonzero resolved tendency.");
            model=WVModel(wvt,shouldUseLinearDynamics=family=="barotropic" || family=="constant-hydrostatic" || family=="constant-nonhydrostatic");
            path=fullfile(testCase.folder,configuration+"-source.nc"); file=model.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true); file.outputTimesForIntegrationPeriod(37,37); file.writeTimeStepToOutputFile(37); model.closeNetCDFFile();
        end
        function [x,y,z]=positionCoordinates(~,wvt)
            dx=wvt.Lx/wvt.Nx; dy=wvt.Ly/wvt.Ny;
            x=[.35*dx;wvt.Lx+.35*dx;2*dx;-.65*dx;wvt.Lx-1e-5*dx;0;wvt.Lx];
            y=[.6*dy;-wvt.Ly+.6*dy;3*dy;wvt.Ly-.4*dy;.75*dy;0;wvt.Ly];
            z=zeros(size(x));
            if ismember('z',wvt.spatialDimensionNames), z=[-wvt.Lz*[.23;.61;.84];wvt.z(1);wvt.z(end);wvt.z(1)-.01*wvt.Lz;wvt.z(end)+.01*wvt.Lz]; end
        end
        function tolerance=samplingTolerance(~,wvt,name,expected)
            scale=max(abs(expected(:)));
            tolerance=2e-11*scale+1e-16;
            if endsWith(name,"_portable_catalog_forcing"), tolerance=2e-11*max(abs(wvt.(name)),[],"all")+realmin; end
            % Retain the independently qualified density-recovery bounds.
            % Interpolation itself is checked separately against the C++ grid.
            if name=="eta_true", tolerance=1e-6*wvt.Lz; end
            if name=="ape", tolerance=1e-5*max(abs(wvt.ape),[],"all")+1e-12; end
            if name=="apv", tolerance=1e-6*max(abs(wvt.apv),[],"all")+1e-12; end
        end
        function writeEvidence(~,suffix,rows)
            report=string(getenv("WV_FIELD_SAMPLING_REPORT"));
            if report~="", writeJSON(replace(report,".json","-"+suffix+".json"),struct(schema="wvm-field-sampling-parity-v1",rows=[rows{:}])); end
        end
        function values=interpolateFullField(~,wvt,x,y,z,method,field)
            dimensions={'x','y'};
            if ~ismatrix(field), dimensions{end+1}='z'; end
            name='portableSamplingGrid';
            annotation=WVVariableAnnotation(name,dimensions,'1','independent full-grid interpolation control');
            operation=WVOperation(name,annotation,@(~)field);
            wvt.addOperation(operation,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            values=wvt.variableAtPositionWithName(x,y,z,name,interpolationMethod=char(method));
        end
        function values=profileColumn(~,wvt,name,xIndex,yIndex)
            value=wvt.(name);
            values=squeeze(value(xIndex,yIndex,:));
            values=values(:);
        end
    end
end
function writeJSON(path,value)
file=fopen(path,"w"); cleanup=onCleanup(@()fclose(file)); fwrite(file,jsonencode(value,PrettyPrint=true));
end
function [status,output]=cleanSystem(command)
[status,output]=system("env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_FALLBACK_LIBRARY_PATH "+command);
end
function value=shellQuote(value)
value="'"+replace(string(value),"'","'""'""'")+"'";
end
