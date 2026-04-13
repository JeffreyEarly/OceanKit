classdef GlobalConstraintUnitTests < matlab.unittest.TestCase

    methods (Test)
        function positiveConstraintHasNoAssociatedDimension(testCase)
            constraint = GlobalConstraint.positive();

            testCase.verifyEqual(constraint.shape, "positive")
            testCase.verifyEmpty(constraint.dimension)
        end

        function monotonicIncreasingDefaultsToFirstDimension(testCase)
            constraint = GlobalConstraint.monotonicIncreasing();

            testCase.verifyEqual(constraint.shape, "monotonicIncreasing")
            testCase.verifyEqual(constraint.dimension, 1)
        end

        function monotonicDecreasingStoresRequestedDimension(testCase)
            constraint = GlobalConstraint.monotonicDecreasing(dimension=2);

            testCase.verifyEqual(constraint.shape, "monotonicDecreasing")
            testCase.verifyEqual(constraint.dimension, 2)
        end

    end
end
