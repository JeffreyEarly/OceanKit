classdef TestSpectralOutputRestart < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function spectralGroupsSurviveInitialAndPopulatedRestarts(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for family = ["constant-hydrostatic","constant-nonhydrostatic","barotropic","stratified-qg","hydrostatic","boussinesq"]
                for multipleGroups = [false true]
                    wvt = TestSpectralOutputRestart.transform(family);
                    model = WVModel(wvt,shouldUseLinearDynamics=true);
                    names = {'Apt','Amt','A0t','u','v','totalEnergy'};
                    names = names(ismember(names,wvt.variableNames));
                    model.eulerianObservingSystem.addNetCDFOutputVariables(names{:});
                    source = fullfile(fixture.Folder,"source.nc");
                    file = model.createNetCDFFileForModelOutput(source,outputInterval=.5,shouldOverwriteExisting=true);
                    groups = "wave-vortex";
                    if multipleGroups
                        dense = file.addNewEvenlySpacedOutputGroup("dense",outputInterval=.125,initialTime=17,finalTime=18);
                        dense.addObservingSystem(WVEulerianFields(model,fieldNames=names));
                        groups(end+1) = "dense";
                    end
                    file.outputTimesForIntegrationPeriod(17,18);
                    file.writeTimeStepToOutputFile(17);
                    model.closeNetCDFFile();
                    paths = fullfile(fixture.Folder,["control.nc","segmented.nc"]);
                    for scenario = 1:2
                        copyfile(source,paths(scenario));
                        times = 18;
                        if scenario == 2, times = [17.5 18]; end
                        for finalTime = times
                            restored = testCase.verifyWarningFree(@()WVModel.modelFromFile(char(paths(scenario))));
                            cleanup = onCleanup(@()restored.closeNetCDFFile());
                            testCase.verifyClass(restored.wvt,class(wvt));
                            for coordinate = ["x","y","k","l","j"]
                                testCase.verifyEqual(restored.wvt.(coordinate),wvt.(coordinate),family+" "+coordinate);
                            end
                            if family ~= "barotropic", testCase.verifyEqual(restored.wvt.z,wvt.z); end
                            testCase.verifyEqual(reshape(restored.outputFiles.outputGroupNames,1,[]),groups);
                            restored.setupIntegrator(integratorType="fixed",deltaT=.25);
                            restored.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                            restored.closeNetCDFFile();
                            clear cleanup
                        end
                    end
                    for group = groups
                        prefix = "/"+group+"/";
                        expectedTimes = (17:(.5-.375*double(group=="dense")):18)';
                        testCase.verifyEqual(ncread(paths(1),prefix+"t"),expectedTimes);
                        testCase.verifyEqual(ncread(paths(2),prefix+"t"),expectedTimes);
                        info = ncinfo(paths(1),"/"+group);
                        for variable = reshape(info.Variables,1,[])
                            if ~ismember(erase(string(variable.Name),["_real","_imag"]),[string(names),"Ap","Am","A0","t","j","kl","x","y","z"]), continue; end
                            name = prefix+string(variable.Name);
                            expected = ncread(paths(1),name);
                            actual = ncread(paths(2),name);
                            testCase.verifyTrue(all(isfinite(actual),"all"));
                            testCase.verifyEqual(actual,expected,family+" "+name,AbsTol=1e-12*max(max(abs(expected),[],"all"),realmin));
                            testCase.verifyEqual(rmfield(ncinfo(paths(2),name),'Filename'),rmfield(ncinfo(paths(1),name),'Filename'));
                        end
                    end
                end
            end
        end

        function unusedOutputCoordinatesCannotChangeConstantGeometry(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            wvt = TestSpectralOutputRestart.transform("constant-hydrostatic");
            path = fullfile(fixture.Folder,"owned-inputs.nc");
            properties = WVGeometryDoublyPeriodicStratifiedConstant.namesOfRequiredPropertiesForGeometry();
            file = wvt.writeToFile(path,properties{:},shouldAddRequiredProperties=false);
            cleanup = onCleanup(@()file.close());
            testCase.verifyFalse(any(string({file.realVariables.name})=="j"));
            for name = ["first","second"]
                output = file.addGroup(name);
                output.addDimension('j',[81;93]);
                output.addVariable('N0',{},99);
            end
            restored = testCase.verifyWarningFree(@()WVTransformConstantStratification.transformFromGroup(file));
            for name = ["j","z","x","N0","Lz","Nz"]
                testCase.verifyEqual(restored.(name),wvt.(name));
            end
            % The generic lookup must still reject ambiguous descendant names.
            testCase.verifyError(@()file.variableWithName('j'),'');
        end
    end

    methods (Static,Access=private)
        function wvt = transform(family)
            switch family
                case "barotropic"
                    wvt = WVTransformBarotropicQG([17000 11000],[8 6],j=1,shouldAntialias=false);
                case "stratified-qg"
                    wvt = WVTransformStratifiedQG([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
                case "boussinesq"
                    wvt = WVTransformBoussinesq([17000 11000 1000],[8 6 9],Nj=4,N2Function=@(z)1e-4*exp(z/700),shouldAntialias=false);
                otherwise
                    wvt = WVTransformConstantStratification([17000 11000 1000],[8 6 9],N0=5.2e-3,isHydrostatic=family=="constant-hydrostatic",shouldAntialias=false);
            end
            wvt.removeAllForcing();
            n = reshape(1:numel(wvt.A0),size(wvt.A0));
            wvt.A0 = 1e-7*complex(sin(.7*n),cos(.3*n)).*(wvt.Kh>0);
            if ~ismember(family,["barotropic","stratified-qg"])
                wvt.Ap = 1e-4*complex(sin(.3*n),cos(.7*n)).*(wvt.Kh>0 & wvt.J>0);
                wvt.Am = 2*wvt.Ap;
            end
            wvt.t = 17;
        end
    end
end
