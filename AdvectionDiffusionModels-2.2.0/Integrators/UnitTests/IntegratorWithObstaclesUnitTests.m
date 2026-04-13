classdef IntegratorWithObstaclesUnitTests < matlab.unittest.TestCase

    methods (Test)
        function reflectsParticleOffObstacleFace(testCase)
            obstacle = polyshape([0 1 1 0],[0 0 1 1]);
            y0 = [-0.5 0.5; -0.5 1.5];
            f = @(~, y) repmat([2 0], size(y,1), 1);

            integrator = IntegratorWithObstacles(f, y0, dt=1, obstacles=obstacle);

            y = integrator.advanceOneStep();

            reflectedExpected = [-1.5 0.5];
            controlExpected = [1.5 1.5];

            testCase.verifyEqual(y(1,:), reflectedExpected, "AbsTol", 2e-3);
            testCase.verifyEqual(y(2,:), controlExpected, "AbsTol", 10*eps);
            testCase.verifyFalse(any(isinterior(obstacle, y(:,1), y(:,2))));
            testCase.verifyEqual(integrator.currentTime, 1);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function handlesParticlesStartingOnObstacleBoundary(testCase)
            obstacle = polyshape([0 1 1 0],[0 0 1 1]);
            y0 = [0 0.25; 0 0.75];
            f = @(~, y) [1 - 2*(y(:,2) >= 0.5), zeros(size(y,1),1)];

            integrator = IntegratorWithObstacles(f, y0, dt=0.5, obstacles=obstacle);

            y = integrator.advanceOneStep();

            reflectedExpected = [-0.5 0.25];
            outwardExpected = [-0.5 0.75];

            testCase.verifyEqual(y(1,:), reflectedExpected, "AbsTol", 2e-3);
            testCase.verifyEqual(y(2,:), outwardExpected, "AbsTol", 10*eps);
            testCase.verifyFalse(any(isinterior(obstacle, y(:,1), y(:,2))));
            testCase.verifyEqual(integrator.currentTime, 0.5);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function handlesParticleTouchingObstacleVertex(testCase)
            obstacle = polyshape([0 1 1 0],[0 0 1 1]);
            f = @(~, y) repmat([1 1], size(y,1), 1);

            integrator = IntegratorWithObstacles(f, [-0.5 -0.5], dt=0.5, obstacles=obstacle);

            y = integrator.advanceOneStep();

            testCase.verifyEqual(y(1), -1e-3, "AbsTol", 2e-3);
            testCase.verifyEqual(y(2), 0, "AbsTol", 1e-12);
            testCase.verifyFalse(isinterior(obstacle, y(1), y(2)));
            testCase.verifyEqual(integrator.currentTime, 0.5);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function reflectsWrappedParticleAgainstObstacle(testCase)
            obstacle = polyshape([0.5 1.5 1.5 0.5],[0 0 1 1]);
            y0 = [3.75 0.5; 3.75 1.5];
            f = @(~, y) repmat([1 0], size(y,1), 1);
            ymin = [0 -Inf];
            ymax = [4 Inf];
            isPeriodic = [true false];

            integrator = IntegratorWithObstacles(f, y0, dt=1, ymin=ymin, ymax=ymax, obstacles=obstacle, isPeriodic=isPeriodic);

            y = integrator.advanceOneStep();
            wrappedX = mod(y(:,1),4);

            testCase.verifyEqual(wrappedX(1), 0.25, "AbsTol", 3e-3);
            testCase.verifyEqual(y(1,2), 0.5, "AbsTol", 10*eps);
            testCase.verifyGreaterThan(y(1,1), 4);
            testCase.verifyEqual(y(2,:), [4.75 1.5], "AbsTol", 10*eps);
            testCase.verifyFalse(any(isinterior(obstacle, wrappedX, y(:,2))));
            testCase.verifyEqual(integrator.currentTime, 1);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function reflectsAcrossMultipleObstaclesInOneStep(testCase)
            obstacles = [polyshape([0 1 1 0],[0 0 1 1]); polyshape([-4 -3 -3 -4],[0 0 1 1])];
            f = @(~, y) repmat([5 0], size(y,1), 1);

            integrator = IntegratorWithObstacles(f, [-0.5 0.5], dt=1, obstacles=obstacles);

            y = integrator.advanceOneStep();

            testCase.verifyEqual(y, [-1.5 0.5], "AbsTol", 2e-3);
            testCase.verifyFalse(any(isinterior(obstacles(1), y(1), y(2))));
            testCase.verifyFalse(any(isinterior(obstacles(2), y(1), y(2))));
            testCase.verifyEqual(integrator.currentTime, 1);
            testCase.verifyEqual(integrator.totalIterations, 1);
        end

        function treatsNearBoundaryStartsContinuously(testCase)
            obstacle = polyshape([0 1 1 0],[0 0 1 1]);
            y0 = [0 0.25; -1e-14 0.25; -1e-12 0.25; -1e-10 0.25];
            f = @(~, y) repmat([1 0], size(y,1), 1);

            integrator = IntegratorWithObstacles(f, y0, dt=0.5, obstacles=obstacle);

            y = integrator.advanceOneStep();

            testCase.verifyEqual(y, repmat([-0.5 0.25], size(y0,1), 1), "AbsTol", 1e-8);
            testCase.verifyFalse(any(isinterior(obstacle, y(:,1), y(:,2))));
        end

        function rejectsOverlappingObstaclesAtConstruction(testCase)
            obstacles = [polyshape([0 1 1 0],[0 0 1 1]); polyshape([0.5 1.5 1.5 0.5],[0 0 1 1])];
            f = @(~, y) zeros(size(y));

            didThrow = false;
            try
                IntegratorWithObstacles(f, [0 0], dt=1, obstacles=obstacles);
            catch exception
                didThrow = true;
                testCase.verifySubstring(exception.message, "overlapping polygons");
            end

            testCase.verifyTrue(didThrow);
        end

        function rejectsInitialParticlesInsideObstacles(testCase)
            obstacle = polyshape([0 1 1 0],[0 0 1 1]);
            f = @(~, y) zeros(size(y));

            didThrow = false;
            try
                IntegratorWithObstacles(f, [1e-10 0.25], dt=0.5, obstacles=obstacle);
            catch exception
                didThrow = true;
                testCase.verifySubstring(exception.message, "inside obstacle");
            end

            testCase.verifyTrue(didThrow);
        end
    end

end
