classdef IntegratorWithObstaclesKnownIssueTests < matlab.unittest.TestCase

    methods (Test)
        function revealsPeriodicMultiWrapObstacleMiss(testCase)
            obstacle = polyshape([0.5 1.5 1.5 0.5],[0 0 1 1]);
            f = @(~, y) repmat([1 0], size(y,1), 1);
            ymin = [0 -Inf];
            ymax = [4 Inf];
            isPeriodic = [true false];

            integrator = IntegratorWithObstacles(f, [2 0.5], dt=9, ymin=ymin, ymax=ymax, obstacles=obstacle, isPeriodic=isPeriodic);

            y = integrator.advanceOneStep();

            % Unwrapped motion should reflect off the obstacle copies at
            % x = 4.5 and x = 1.5, and then off x = 4.5 again, ending at x = 4.
            testCase.verifyEqual(y, [4 0.5], "AbsTol", 2e-3);
            testCase.verifyFalse(isinterior(obstacle, mod(y(1),4), y(2)));
            testCase.verifyEqual(integrator.currentTime, 9);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function revealsReflectionLoopBailout(testCase)
            obstacles = [polyshape([0 1 1 0],[0 0 1 1]); polyshape([2 3 3 2],[0 0 1 1])];
            f = @(~, y) repmat([100 0], size(y,1), 1);

            integrator = IntegratorWithObstacles(f, [1.5 0.5], dt=1, obstacles=obstacles);

            lastwarn('');
            y = integrator.advanceOneStep();
            warningMessage = lastwarn;

            testCase.verifyEmpty(warningMessage);
            testCase.verifyEqual(y, [1.5 0.5], "AbsTol", 2e-3);
            testCase.verifyFalse(isinterior(obstacles(1), y(1), y(2)));
            testCase.verifyFalse(isinterior(obstacles(2), y(1), y(2)));
            testCase.verifyEqual(integrator.currentTime, 1);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end
    end

end
