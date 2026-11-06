classdef IMZeroAPVEnergyMatrixTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function indefiniteAndRepeatedEnergyHaveRealOrthonormalModes(testCase)
            for H={diag([2 -3]),2*eye(2),zeros(2),[2 1;1 -2]}
                matrix=H{1};
                [C,depth,evidence]=IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(matrix,.5,"IMEnergyAudit");
                testCase.verifyEqual(C.'*C,eye(2),AbsTol=3e-14);
                testCase.verifyEqual(matrix*C,C*diag(depth/(2*.5^2)),AbsTol=3e-14);
                testCase.verifyLessThan(evidence.endpointOrthogonalityError,3e-14);
            end
        end
        function physicalUnitsDoNotChangeTheOrthogonalDirections(testCase)
            H=[3 .2;.2 -1];
            [reference,depth]=IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(H,1,"IMEnergyAudit");
            for scale=[1e-250 1 1e250]
                [C,current]=IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(scale*H,1,"IMEnergyAudit");
                testCase.verifyEqual(C,reference,AbsTol=3e-14);
                testCase.verifyEqual(current/scale,depth,RelTol=3e-13);
            end
        end
        function emptyPagesHaveAnEmptyGeometry(testCase)
            [C,depth,evidence]=IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(zeros(0,0,2),[1 2],"IMEnergyAudit");
            testCase.verifySize(C,[0 0 2]);
            testCase.verifySize(depth,[0 2]);
            testCase.verifyEqual(evidence.endpointOrthogonalityError,[0 0],AbsTol=0);
        end
        function asymmetricFormsFailInTheirOwnUnits(testCase)
            testCase.verifyError(@() IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(1e-250*[1 1;0 1],1,"IMEnergyAudit"),"IMEnergyAudit:NonSymmetricGeneralizedEnergy");
        end
    end
end
