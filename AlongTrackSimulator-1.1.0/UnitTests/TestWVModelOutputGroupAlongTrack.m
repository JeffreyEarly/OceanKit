classdef TestWVModelOutputGroupAlongTrack < matlab.unittest.TestCase

    properties
        wvt
        model
        ats
        fieldNames
    end

    methods (TestClassSetup)
        function createModel(testCase)
            Lxyz = [2000e3 1000e3 500];
            Nxyz = [8 8 5];
            N2 = @(z) (5.2e-3)^2*ones(size(z));
            testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=24,shouldAntialias=false);

            rng(1);
            testCase.wvt.initWithRandomFlow(uvMax=0.01);
            testCase.wvt.t = 3600;

            balancedComponent = testCase.wvt.geostrophicComponent + testCase.wvt.mdaComponent;
            balancedComponent.name = "balanced";
            balancedComponent.shortName = "balanced";
            balancedComponent.abbreviatedName = "balanced";
            testCase.wvt.addFlowComponent(balancedComponent);

            sshAnnotation = testCase.wvt.propertyAnnotationWithName("ssh");
            sshAnnotation.attributes('test_attribute') = "preserved";
            testCase.wvt.addForcing(WVAdaptiveDamping(testCase.wvt));
            testCase.model = WVModel(testCase.wvt);
            testCase.ats = AlongTrackSimulator();
            testCase.fieldNames = ["ssh","ssu","ssv", ...
                "ssh_balanced","ssu_balanced","ssv_balanced", ...
                "ssh_w","ssu_w","ssv_w", ...
                "ssh_io","ssu_io","ssv_io"];
        end
    end

    methods (Test)
        function testDefaultOutputIsSSHOnly(testCase)
            group = testCase.outputGroupWithTestTrack();
            testCase.verifyEqual(group.fieldNames,"ssh");

            [ncfile,fileCleanup] = testCase.temporaryNetCDFFile(); %#ok<ASGLU>
            group.writeTimeStepToNetCDFFile(ncfile,testCase.wvt.t);
            missionGroup = ncfile.groupWithName("alg");

            testCase.verifyTrue(missionGroup.hasVariableWithName("ssh"));
            testCase.verifyFalse(missionGroup.hasVariableWithName("ssu"));
            testCase.verifyFalse(missionGroup.hasVariableWithName("ssv"));
        end

        function testRequestedFieldsAreWrittenAlongTrack(testCase)
            group = testCase.outputGroupWithTestTrack(fieldNames=testCase.fieldNames);
            track = group.tracks{1};
            expectedValues = cell(size(testCase.fieldNames));
            fieldNamesCell = cellstr(testCase.fieldNames);
            [expectedValues{:}] = testCase.wvt.variableAtPositionWithName(track.x,track.y,[],fieldNamesCell{:});

            [ncfile,fileCleanup] = testCase.temporaryNetCDFFile(); %#ok<ASGLU>
            group.writeTimeStepToNetCDFFile(ncfile,testCase.wvt.t);
            missionGroup = ncfile.groupWithName("alg");

            testCase.verifyEqual(missionGroup.readVariables("track_x"),reshape(track.x,[],1),AbsTol=1e-12);
            testCase.verifyEqual(missionGroup.readVariables("track_y"),reshape(track.y,[],1),AbsTol=1e-12);
            testCase.verifyEqual(missionGroup.readVariables("t"),reshape(track.t,[],1),AbsTol=1e-12);

            for iField = 1:length(testCase.fieldNames)
                fieldName = testCase.fieldNames(iField);
                testCase.verifyTrue(missionGroup.hasVariableWithName(char(fieldName)));
                actual = missionGroup.readVariables(char(fieldName));
                testCase.verifyEqual(actual,reshape(expectedValues{iField},[],1),AbsTol=1e-12);

                annotation = testCase.wvt.propertyAnnotationWithName(fieldName);
                variable = missionGroup.variableWithName(char(fieldName));
                testCase.verifyEqual(string(variable.attributes('units')),string(annotation.units));
                testCase.verifyTrue(contains(string(variable.attributes('long_name')),string(annotation.description)));
                testCase.verifyTrue(contains(string(variable.attributes('long_name')),"sampled along the alg mission ground track"));
            end
            testCase.verifyEqual(string(missionGroup.variableWithName("ssh").attributes('test_attribute')),"preserved");
        end

        function testConvenienceMethodAddsMissionGroups(testCase)
            outputPath = string(tempname) + ".nc";
            outputFile = WVModelOutputFile(testCase.model,outputPath);
            groups = testCase.ats.addMissionsToOutputFile(outputFile,missionNames="alg",fieldNames=testCase.fieldNames);

            testCase.verifySize(groups,[1 1]);
            testCase.verifyEqual(groups.fieldNames,testCase.fieldNames);
            testCase.verifyEqual(outputFile.outputGroupNames,"alg");
            testCase.verifyEqual(outputFile.outputGroupWithName('alg'),groups);
        end

        function testRejectsInvalidFieldSelections(testCase)
            constructor = @(fieldNames) WVModelOutputGroupAlongTrack(testCase.model,"alg",testCase.ats,fieldNames=fieldNames);
            testCase.verifyError(@() constructor(strings(1,0)),"WVModelOutputGroupAlongTrack:EmptyFieldNames");
            testCase.verifyError(@() constructor(["ssh","ssh"]),"WVModelOutputGroupAlongTrack:DuplicateFieldNames");
            testCase.verifyError(@() constructor("notARegisteredField"),"WVModelOutputGroupAlongTrack:UnknownField");
            testCase.verifyError(@() constructor("u"),"WVModelOutputGroupAlongTrack:InvalidFieldDimensions");
        end

        function testLegacyRepeatMissionHelper(testCase)
            group = AlongTrackSimulator.wvmOutputGroupForRepeatMissionWithName(testCase.model,"s6a");

            testCase.verifyClass(group,"WVModelOutputGroupAlongTrack");
            testCase.verifyEqual(group.missionName,"s6a");
            testCase.verifyEqual(group.fieldNames,"ssh");
            testCase.verifyNotEmpty(group.tracks);
        end
    end

    methods (Access=private)
        function group = outputGroupWithTestTrack(testCase,options)
            arguments
                testCase
                options.fieldNames string = "ssh"
            end
            group = WVModelOutputGroupAlongTrack(testCase.model,"alg",testCase.ats,fieldNames=options.fieldNames);
            dx = testCase.wvt.x(2)-testCase.wvt.x(1);
            dy = testCase.wvt.y(2)-testCase.wvt.y(1);
            track = struct;
            track.x = [3.25*dx -0.25*dx 3.25*dx testCase.wvt.Lx+0.2*dx];
            track.y = [3.4*dy 3.4*dy testCase.wvt.Ly+0.3*dy testCase.wvt.Ly-0.25*dy];
            track.t = testCase.wvt.t + (0:3);
            group.tracks = {track};
            group.firstPassoverTime = testCase.wvt.t;
        end

        function [ncfile,fileCleanup] = temporaryNetCDFFile(~)
            outputPath = string(tempname) + ".nc";
            ncfile = NetCDFFile(outputPath);
            fileCleanup = onCleanup(@() TestWVModelOutputGroupAlongTrack.closeAndDelete(ncfile,outputPath));
        end
    end

    methods (Static,Access=private)
        function closeAndDelete(ncfile,outputPath)
            ncfile.close();
            if isfile(outputPath)
                delete(outputPath);
            end
        end
    end

end
