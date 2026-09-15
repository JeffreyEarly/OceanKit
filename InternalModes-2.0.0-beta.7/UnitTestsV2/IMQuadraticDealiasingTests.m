classdef IMQuadraticDealiasingTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addPackage(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
        end
    end

    methods (Test)
        function policiesOperateOnTheFullLinearInventory(testCase)
            F = eye(10);
            fixed = assessQuadraticDealiasing(F,[],9,quadraticDealiasing="fixedFraction");
            testCase.verifyEqual(fixed.accepted,[true(6,1);false(4,1)]);
            testCase.verifyTrue(all(isnan(fixed.effectiveDegree)));
            none = assessQuadraticDealiasing(F,[],9,quadraticDealiasing="none");
            testCase.verifyEqual(none.accepted,true(10,1));
            testCase.verifyTrue(all(isnan(none.tailEnergyFractionF)));
            empty = assessQuadraticDealiasing(F,[],9,quadraticDealiasing="fixedFraction",retainedFraction=0);
            full = assessQuadraticDealiasing(F,[],9,quadraticDealiasing="fixedFraction",retainedFraction=1);
            testCase.verifyFalse(any(empty.accepted));
            testCase.verifyTrue(all(full.accepted));
        end

        function pureDegreesGiveIndependentModeMasks(testCase)
            coefficients = eye(13);
            order = [1 10 3 9 13];
            result = assessQuadraticDealiasing(coefficients(:,order),[],12,representation="coefficients",quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(result.effectiveDegree,[0;9;2;8;12]);
            testCase.verifyEqual(result.accepted,[true;false;true;true;false]);
            testCase.verifyEqual(result.bandwidthLimit,8);
            testCase.verifyEqual(result.tailEnergyFractionF,[0;1;0;0;1]);
        end

        function chebyshevWeightedEnergySetsTheCutoff(testCase)
            coefficients = [1;zeros(4,1);1];
            low = assessQuadraticDealiasing(coefficients,[],4,representation="coefficients",energyFraction=2/3,quadraticDealiasing="effectiveBandwidth");
            high = assessQuadraticDealiasing(coefficients,[],4,representation="coefficients",energyFraction=.9,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(low.effectiveDegree,0);
            testCase.verifyEqual(high.effectiveDegree,5);
            testCase.verifyEqual(high.tailEnergyFractionF,1/3,AbsTol=eps);
            complete = assessQuadraticDealiasing(coefficients,[],4,representation="coefficients",energyFraction=1,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(complete.effectiveDegree,5);
        end

        function eitherPhysicalChannelCanLimitAMode(testCase)
            F = zeros(17,3);
            G = F;
            F(1,1) = 1;
            F(3,2) = 1;
            G(15,2) = 1;
            G(5,3) = 1;
            result = assessQuadraticDealiasing(F,G,12,representation="coefficients",quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(result.effectiveDegreeF,[0;2;0]);
            testCase.verifyEqual(result.effectiveDegreeG,[0;14;4]);
            testCase.verifyEqual(result.effectiveDegree,[0;14;4]);
            testCase.verifyEqual(result.accepted,[true;false;true]);
        end

        function fullCoveragePreservesTinyNonzeroTails(testCase)
            coefficients = [1 1 1e300;1e-10 1e-200 1e-300];
            result = assessQuadraticDealiasing(coefficients,[],0,representation="coefficients",energyFraction=1,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(result.effectiveDegree,ones(3,1));
            testCase.verifyFalse(any(result.accepted));
            testCase.verifyEqual(result.tailEnergyFractionF(1),5e-21,RelTol=4*eps);
            testCase.verifyEqual(result.tailEnergyFractionF(2:3),zeros(2,1));
            relaxed = assessQuadraticDealiasing(coefficients(:,1),[],0,representation="coefficients",energyFraction=.99,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(relaxed.effectiveDegree,0);
            testCase.verifyTrue(relaxed.accepted);
            testCase.verifyEqual(relaxed.tailEnergyFractionF,5e-21,RelTol=4*eps);
        end

        function scalingIsIndependentAcrossModesAndChannels(testCase)
            xi = -cos(pi*(0:64)'/64);
            F = [ones(size(xi)),cos(6*acos(xi)),cos(15*acos(xi))];
            G = [cos(3*acos(xi)),cos(7*acos(xi)),cos(2*acos(xi))];
            expected = assessQuadraticDealiasing(F,G,18,quadraticDealiasing="effectiveBandwidth");
            actual = assessQuadraticDealiasing(F.*[1e-280 -1e280 3],G.*[-1e260 1e-260 -2],18,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(actual.effectiveDegree,expected.effectiveDegree);
            testCase.verifyEqual(actual.accepted,expected.accepted);
            testCase.verifyEqual(actual.tailEnergyFractionF,expected.tailEnergyFractionF,AbsTol=2e-15);
            testCase.verifyEqual(actual.tailEnergyFractionG,expected.tailEnergyFractionG,AbsTol=2e-15);
            tiny = assessQuadraticDealiasing(realmin*eps,[],0,representation="coefficients",quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(tiny.effectiveDegree,0);
            testCase.verifyTrue(tiny.accepted);
        end

        function valueTransformAgreesWithKnownCoefficients(testCase)
            Q = 65;
            xi = -cos(pi*(0:Q-1)'/(Q-1));
            coefficients = zeros(Q,3);
            coefficients([1 6 18],1) = [1 .2 .5];
            coefficients([2 9 32],2) = [-.5 1 .3];
            coefficients([1 65],3) = [1 2];
            values = cos(acos(xi)*(0:Q-1))*coefficients;
            direct = assessQuadraticDealiasing(coefficients,[],32,representation="coefficients",quadraticDealiasing="effectiveBandwidth");
            sampled = assessQuadraticDealiasing(values,[],32,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(sampled.effectiveDegree,direct.effectiveDegree);
            testCase.verifyEqual(sampled.accepted,direct.accepted);
            testCase.verifyEqual(sampled.tailEnergyFractionF,direct.tailEnergyFractionF,AbsTol=3e-14);
            reflected = assessQuadraticDealiasing(flipud(values),[],32,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(reflected.effectiveDegree,direct.effectiveDegree);
            testCase.verifyEqual(reflected.tailEnergyFractionF,direct.tailEnergyFractionF,AbsTol=3e-14);
        end

        function nativeCoordinatesMustBeEvaluatedOnTheCommonGrid(testCase)
            xi = -cos(pi*(0:64)'/64);
            z = (xi+.2*xi.^3)/1.2;
            fromPhysicalCoordinate = 2*z.^2-1;
            commonCoefficients = [-23/288;0;575/576;0;23/288;0;1/576];
            fromCommonCoordinate = cos(acos(xi)*(0:6))*commonCoefficients;
            physical = assessQuadraticDealiasing(fromPhysicalCoordinate,[],8,energyFraction=.999999,quadraticDealiasing="effectiveBandwidth");
            common = assessQuadraticDealiasing(fromCommonCoordinate,[],8,energyFraction=.999999,quadraticDealiasing="effectiveBandwidth");
            native = assessQuadraticDealiasing([0;0;1],[],8,representation="coefficients",energyFraction=.999999,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(physical.effectiveDegree,common.effectiveDegree);
            testCase.verifyEqual(physical.tailEnergyFractionF,common.tailEnergyFractionF,AbsTol=2e-14);
            testCase.verifyEqual(common.effectiveDegree,6);
            testCase.verifyEqual(native.effectiveDegree,2);
            testCase.verifyFalse(common.accepted);
            testCase.verifyTrue(native.accepted);
        end

        function batchingPreservesMeasuredResults(testCase)
            Q = 33;
            xi = -cos(pi*(0:Q-1)'/(Q-1));
            F = cos(acos(xi)*(0:25));
            G = sin(acos(xi)*(0:25));
            scalar = assessQuadraticDealiasing(F,G,24,batchSize=1,quadraticDealiasing="effectiveBandwidth");
            uneven = assessQuadraticDealiasing(F,G,24,batchSize=7,quadraticDealiasing="effectiveBandwidth");
            together = assessQuadraticDealiasing(F,G,24,batchSize=128,quadraticDealiasing="effectiveBandwidth");
            for actual = [uneven,together]
                testCase.verifyEqual(actual.effectiveDegreeF,scalar.effectiveDegreeF);
                testCase.verifyEqual(actual.effectiveDegreeG,scalar.effectiveDegreeG);
                testCase.verifyEqual(actual.accepted,scalar.accepted);
                testCase.verifyEqual(actual.tailEnergyFractionF,scalar.tailEnergyFractionF,AbsTol=2e-14);
                testCase.verifyEqual(actual.tailEnergyFractionG,scalar.tailEnergyFractionG,AbsTol=2e-14);
            end
        end

        function absentNullEmptyAndConstantInventoriesAreDefined(testCase)
            for policy = ["none","fixedFraction","effectiveBandwidth"]
                empty = assessQuadraticDealiasing([],[],0,quadraticDealiasing=policy);
                testCase.verifySize(empty.accepted,[0 1]);
                sizedEmpty = assessQuadraticDealiasing(zeros(17,0),[],8,quadraticDealiasing=policy);
                testCase.verifySize(sizedEmpty.effectiveDegree,[0 1]);
                testCase.verifyEqual(sizedEmpty.sampleCount,17);
            end
            zerosOnly = assessQuadraticDealiasing([],zeros(17,3),0,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(zerosOnly.effectiveDegree,zeros(3,1));
            testCase.verifyEqual(zerosOnly.tailEnergyFractionG,zeros(3,1));
            testCase.verifyTrue(all(zerosOnly.accepted));
            constants = assessQuadraticDealiasing([0 2 -3],[],0,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(constants.effectiveDegree,zeros(3,1));
            twoPoints = assessQuadraticDealiasing([-1;1],[],1,bandwidthFraction=1,quadraticDealiasing="effectiveBandwidth");
            testCase.verifyEqual(twoPoints.effectiveDegree,1);
            testCase.verifyTrue(twoPoints.accepted);
        end

        function malformedShapesFailClearly(testCase)
            testCase.verifyError(@() assessQuadraticDealiasing(zeros(3,2),zeros(4,2),2,quadraticDealiasing="effectiveBandwidth"),"assessQuadraticDealiasing:ChannelSizeMismatch");
            testCase.verifyError(@() assessQuadraticDealiasing(zeros(3,2),zeros(3,1),2,quadraticDealiasing="effectiveBandwidth"),"assessQuadraticDealiasing:ChannelSizeMismatch");
            testCase.verifyError(@() assessQuadraticDealiasing(zeros(0,3),[],2,quadraticDealiasing="effectiveBandwidth"),"assessQuadraticDealiasing:MissingSamples");
        end
    end
end
