classdef TestStratifiedModalRecord < matlab.unittest.TestCase
    properties (SetAccess = private)
        temporaryFolder (1,1) string
        reader (1,1) string
    end

    methods (TestClassSetup)
        function buildReader(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
            testCase.reader = string(getenv("WV_MODAL_DUMP_EXECUTABLE"));
            if testCase.reader == ""
                build = fullfile(testCase.temporaryFolder,"build");
                [status,output] = systemWithoutMatlabRuntime("cmake -S " + shellQuote(fullfile(root,"PortableRuntime")) + " -B " + shellQuote(build) + " -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DWV_ENABLE_ACCELERATE=OFF");
                testCase.assertEqual(status,0,output)
                [status,output] = systemWithoutMatlabRuntime("cmake --build " + shellQuote(build) + " --parallel 4 --target WVStratifiedModalDump");
                testCase.assertEqual(status,0,output)
                testCase.reader = fullfile(build,"WVStratifiedModalDump");
            end
            testCase.assertTrue(isfile(testCase.reader))
        end
    end

    methods (Test, TestTags = "full")
        function olderFilesRemainMatlabCompatible(testCase)
            constructors = {@WVTransformStratifiedQG,@WVTransformHydrostatic};
            for iTransform = 1:numel(constructors)
                wvt = constructors{iTransform}([17000 23000 1000],[9 8 9],Nj=5,N2Function=@(z) 1e-4*exp(z/800),shouldAntialias=false);
                path = fullfile(testCase.temporaryFolder,"older.nc");
                % Reproduce the original default schema from its unchanged
                % required-property list, without the additive writer fields.
                required = wvt.requiredProperties;
                testCase.verifyFalse(any(ismember({'N2','rho_nm0','k','l'},required)))
                file = wvt.writeToFile(char(path),required{:},shouldAddRequiredProperties=false,shouldOverwriteExisting=true);
                file.close();
                info = ncinfo(path);
                testCase.verifyFalse(any(ismember(["N2","rho_nm0","k","l"],string({info.Variables.Name}))))
                [restored,file] = WVTransform.waveVortexTransformFromFile(char(path));
                cleanup = onCleanup(@()file.close());
                geometry = WVGeometryDoublyPeriodicStratified.geometryFromGroup(file);
                testCase.verifyClass(restored,class(wvt))
                for name = ["j","z","N2","rho_nm0","PF0inv","PF0","QG0inv","QG0","P0","Q0","h_0","dLnN2","z_int"]
                    tolerance = 2e-11*max(1,max(abs(wvt.(name)),[],"all"));
                    testCase.verifyEqual(restored.(name),wvt.(name),name,AbsTol=tolerance)
                    testCase.verifyEqual(geometry.(name),wvt.(name),name,AbsTol=tolerance)
                end
                clear cleanup
            end
        end

        function defaultModelOutputContainsScientificRecords(testCase)
            constructors = {@WVTransformStratifiedQG,@WVTransformHydrostatic};
            for iTransform = 1:numel(constructors)
                wvt = constructors{iTransform}([17000 23000 1000],[9 8 9],Nj=5,N2Function=@(z) 1e-4*exp(z/800),shouldAntialias=false);
                model = WVModel(wvt,shouldUseLinearDynamics=true);
                cleanup = onCleanup(@()model.closeNetCDFFile());
                path = fullfile(testCase.temporaryFolder,"model.nc");
                model.createNetCDFFileForModelOutput(char(path),outputInterval=1,shouldOverwriteExisting=true);
                model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                model.closeNetCDFFile();
                clear cleanup
                [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.reader) + " " + shellQuote(path));
                testCase.assertEqual(status,0,text)
                report = jsondecode(text);
                for name = ["N2","rho_nm0","k","l"]
                    testCase.verifyEqual(report.(name),wvt.(name))
                end
                restored = WVModel.modelFromFile(char(path));
                cleanup = onCleanup(@()restored.closeNetCDFFile());
                testCase.verifyClass(restored.wvt,class(wvt))
                testCase.verifyEqual(restored.t,1)
                clear cleanup
            end
        end

        function explicitPropertyOnlySavesRemainUnchanged(testCase)
            wvt = WVTransformStratifiedQG([17000 23000 1000],[9 8 9],Nj=5,N2Function=@(z) 1e-4*exp(z/800));
            path = fullfile(testCase.temporaryFolder,"properties.nc");
            file = wvt.writeToFile(char(path),'N2',shouldAddRequiredProperties=false);
            file.close();
            info = ncinfo(path);
            names = string({info.Variables.Name});
            testCase.verifyTrue(ismember("N2",names))
            testCase.verifyFalse(any(ismember(["rho_nm0","k","l","PF0"],names)))
            testCase.verifyEqual(ncread(path,'N2'),wvt.N2)
        end

        function nonprefixModalSubsetPreservesScientificKeys(testCase)
            wvt = WVTransformStratifiedQG([19000 11000 1000],[9 7 10],Nj=6,N2Function=@(z) 1e-4*exp(z/700));
            keep = [1 3 5];
            path = fullfile(testCase.temporaryFolder,"subset.nc");
            dimensions = {"x","y","z","j","kl"};
            for iDimension = 1:numel(dimensions)
                name = dimensions{iDimension};
                values = wvt.(name);
                if name == "j", values = values(keep); end
                nccreate(path,name,"Dimensions",{char(name),numel(values)},"Format","netcdf4");
                ncwrite(path,name,values);
            end
            for name = ["Lx","Ly","Lz","g","rho0","latitude","rotationRate","planetaryRadius"]
                nccreate(path,name,"Format","netcdf4"); ncwrite(path,name,wvt.(name));
            end
            nccreate(path,"shouldAntialias","Datatype","uint8","Format","netcdf4");
            ncwrite(path,"shouldAntialias",uint8(wvt.shouldAntialias));
            fields = {"k",{"kl",wvt.Nkl}; "l",{"kl",wvt.Nkl}; "N2",{"z",wvt.Nz}; "rho_nm0",{"z",wvt.Nz}; "dLnN2",{"z",wvt.Nz}; "z_int",{"z",wvt.Nz}; "P0",{"j",numel(keep)}; "Q0",{"j",numel(keep)}; "h_0",{"j",numel(keep)}; "PF0inv",{"z",wvt.Nz,"j",numel(keep)}; "QG0inv",{"z",wvt.Nz,"j",numel(keep)}; "PF0",{"j",numel(keep),"z",wvt.Nz}; "QG0",{"j",numel(keep),"z",wvt.Nz}};
            for iField = 1:size(fields,1)
                name = fields{iField,1}; values = wvt.(name);
                if any(name == ["P0","Q0","h_0"]), values = values(keep); end
                if any(name == ["PF0inv","QG0inv"]), values = values(:,keep); end
                if any(name == ["PF0","QG0"]), values = values(keep,:); end
                nccreate(path,name,"Dimensions",fields{iField,2},"Format","netcdf4"); ncwrite(path,name,values);
            end
            ncwriteatt(path,"/","model_version",char(wvt.version));
            ncwriteatt(path,"/","WVTransform",class(wvt));
            ncwriteatt(path,"/","AnnotatedClass",class(wvt));
            [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.reader) + " " + shellQuote(path));
            testCase.assertEqual(status,0,text)
            report = jsondecode(text);
            testCase.verifyEqual(report.j,wvt.j(keep))
            testCase.verifyEqual(reshape(report.PF0inv,wvt.Nz,[]),wvt.PF0inv(:,keep))
            F = reshape(report.PF0,numel(keep),[])./report.P0;
            Finv = reshape(report.PF0inv,wvt.Nz,[]).*report.P0';
            testCase.verifyEqual(F*Finv,eye(numel(keep)),AbsTol=2e-10)
            input = 0.1*(1:numel(keep))'+0.03*(0:wvt.Nkl-1)+1i*(-0.2*(1:numel(keep))'+0.01*(0:wvt.Nkl-1));
            expected = wvt.FinvMatrix(:,keep)*input;
            actual = reshape(report.Finv.interleaved.real+1i*report.Finv.interleaved.imag,wvt.Nz,[]);
            testCase.verifyEqual(actual,expected,AbsTol=2e-12)
        end

        function sharedRecordsRoundTripThroughMatlab(testCase)
            constructors = {@WVTransformStratifiedQG,@WVTransformHydrostatic};
            for iTransform = 1:numel(constructors)
                wvt = constructors{iTransform}([17000 23000 1000],[9 8 9],Nj=5,N2Function=@(z) 1e-4*exp(z/800),shouldAntialias=false);
                path = fullfile(testCase.temporaryFolder,"roundtrip.nc");
                file = wvt.writeToFile(char(path),shouldOverwriteExisting=true);
                file.close();
                [restored,file] = WVTransform.waveVortexTransformFromFile(char(path));
                cleanup = onCleanup(@()file.close());
                testCase.verifyClass(restored,class(wvt))
                for name = ["j","z","N2","rho_nm0","PF0inv","PF0","QG0inv","QG0","P0","Q0","h_0","dLnN2","z_int"]
                    testCase.verifyEqual(restored.(name),wvt.(name),name,AbsTol=2e-11*max(1,max(abs(wvt.(name)),[],"all")))
                end
                clear cleanup
                [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.reader) + " " + shellQuote(path));
                testCase.assertEqual(status,0,text)
                report = jsondecode(text);
                testCase.verifyEqual(report.N2,wvt.N2)
                % State/time edits do not alter scientific modal values.
                ncwrite(path,"t",123);
                ncwrite(path,"t0",17);
                [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.reader) + " " + shellQuote(path));
                testCase.assertEqual(status,0,text)
                testCase.verifyEqual(jsondecode(text),report)
            end
        end

        function scientificRecordsAndPreparedOperatorsMatchMatlab(testCase)
            profiles = {@(z) 1e-4*exp(z/700), @(z) 5e-5*(1+0.4*tanh((z+400)/180))};
            sizes = [8 6 9 4; 9 7 10 6; 8 7 7 1];
            for iProfile = 1:numel(profiles)
                for iSize = 1:size(sizes,1)
                    for shouldAntialias = [false true]
                        sz = sizes(iSize,:);
                        wvt = WVTransformStratifiedQG([19000 11000 1000],sz(1:3),Nj=sz(4),N2Function=profiles{iProfile},shouldAntialias=shouldAntialias);
                        testCase.verifyGreaterThan(max(diff(wvt.z))-min(diff(wvt.z)),1)
                        path = fullfile(testCase.temporaryFolder,"modal.nc");
                        file = wvt.writeToFile(char(path),shouldOverwriteExisting=true);
                        file.close();
                        [status,text] = systemWithoutMatlabRuntime(shellQuote(testCase.reader) + " " + shellQuote(path));
                        testCase.assertEqual(status,0,text)
                        report = jsondecode(text);
                        testCase.verifyEqual(string(report.contract),"wave-vortex-stratified-modal-record-v1")
                        testCase.verifyEqual(report.matrixBytes,4*wvt.Nz*wvt.Nj*8)
                        for name = ["x","y","z","j","k","l","N2","rho_nm0","dLnN2","P0","Q0","h_0","z_int","PF0inv","QG0inv","PF0","QG0"]
                            expected = wvt.(name);
                            testCase.verifyEqual(report.(name)(:),expected(:),name)
                        end
                        operators = {"Finv",wvt.FinvMatrix; "F",wvt.FMatrix; "Ginv",wvt.GinvMatrix; "G",wvt.GMatrix; "GToF",wvt.FMatrix*wvt.GinvMatrix; "FToG",wvt.GMatrix*wvt.FinvMatrix};
                        for iOperator = 1:size(operators,1)
                            name = operators{iOperator,1};
                            matrix = operators{iOperator,2};
                            rows = (1:size(matrix,2))';
                            columns = 0:wvt.Nkl-1;
                            input = (0.1*rows+0.03*columns) + 1i*(-0.2*rows+0.01*columns);
                            expected = matrix*input;
                            for representation = ["split","interleaved"]
                                values = report.(name).(representation);
                                actual = reshape(values.real+1i*values.imag,size(expected));
                                testCase.verifyLessThanOrEqual(max(abs(actual-expected),[],"all"),2e-12*max(1,max(abs(expected),[],"all")),name+" "+representation)
                            end
                        end
                        testCase.verifyEqual(wvt.FMatrix*wvt.FinvMatrix,eye(wvt.Nj),AbsTol=2e-10)
                        expectedG = eye(wvt.Nj); expectedG(1,1) = 0;
                        testCase.verifyEqual(wvt.GMatrix*wvt.GinvMatrix,expectedG,AbsTol=2e-10)
                        testCase.verifyEqual(sum(report.z_int),wvt.Lz,RelTol=1e-12)
                        % Reconstruct calculus from decoded scientific records and
                        % compare the current MATLAB operations, including their
                        % barotropic removal and bottom-zero integration constant.
                        PF = reshape(report.PF0,wvt.Nj,wvt.Nz);
                        QG = reshape(report.QG0,wvt.Nj,wvt.Nz);
                        PFinv = reshape(report.PF0inv,wvt.Nz,wvt.Nj);
                        QGinv = reshape(report.QG0inv,wvt.Nz,wvt.Nj);
                        p = report.P0; q = report.Q0; h = report.h_0; N2 = report.N2;
                        DzF = -(N2/wvt.g).*QGinv*((q./p).*PF);
                        DzG = PFinv*((p./(q.*h)).*QG);
                        IntF = QGinv*((h.*q./p).*PF);
                        IntG = -wvt.g*(PFinv-PFinv(1,:))*((p./q).*(QG./N2'));
                        field = reshape(sin(0.17*(1:wvt.Nx*wvt.Ny*wvt.Nz)),wvt.spatialMatrixSize);
                        vertical = reshape(permute(field,[3 1 2]),wvt.Nz,[]);
                        methods = {DzF,wvt.diffZF(field); DzG,wvt.diffZG(field); IntF,wvt.intZF(field); IntG,wvt.intZG(field)};
                        for iMethod = 1:size(methods,1)
                            expected = reshape(permute(methods{iMethod,2},[3 1 2]),wvt.Nz,[]);
                            actual = methods{iMethod,1}*vertical;
                            testCase.verifyEqual(actual,expected,AbsTol=2e-11*max(1,max(abs(expected),[],"all")))
                        end
                        n = 0:wvt.Nx*wvt.Ny*wvt.Nj-1;
                        grid = reshape(sin(0.07*n)+0.2*cos(0.19*n),wvt.Nx,wvt.Ny,wvt.Nj);
                        coefficients = reshape(fft(fft(grid,[],1),[],2),wvt.Nx*wvt.Ny,wvt.Nj)/(wvt.Nx*wvt.Ny);
                        indices = 1+mod(round(wvt.k*wvt.Lx/(2*pi)),wvt.Nx)+wvt.Nx*mod(round(wvt.l*wvt.Ly/(2*pi)),wvt.Ny);
                        expected = coefficients(indices,:).';
                        actual = reshape(report.horizontal.real+1i*report.horizontal.imag,wvt.Nj,wvt.Nkl);
                        testCase.verifyEqual(actual,expected,AbsTol=2e-13)
                    end
                end
            end
        end
    end
end

function value = shellQuote(value)
value = "'" + replace(string(value),"'","'""'""'") + "'";
end

function [status,output] = systemWithoutMatlabRuntime(command)
if isunix && ~ismac
    command = "env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH " + command;
end
[status,output] = system(command);
end
