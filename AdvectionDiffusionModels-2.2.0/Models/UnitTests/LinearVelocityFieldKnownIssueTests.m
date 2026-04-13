classdef LinearVelocityFieldKnownIssueTests < matlab.unittest.TestCase

    methods (Test)
        function particlePathMissesBackgroundVelocityInMatchedFlow(testCase)
            sigma = 8e-6;
            theta = 0;
            zeta = 8e-6;
            u0 = 0.08;
            v0 = -0.05;

            [x0, y0] = LinearVelocityFieldKnownIssueTests.legacyInitialPositions();
            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);

            maxError = LinearVelocityFieldKnownIssueTests.maxPathError(model,x0,y0,2*86400,360,u0,v0);

            testCase.verifyLessThanOrEqual(maxError, 5e-6);
        end

        function particlePathMissesBackgroundVelocityInPureTranslation(testCase)
            u0 = 0.08;
            v0 = -0.05;

            [x0, y0] = LinearVelocityFieldKnownIssueTests.legacyInitialPositions();
            model = LinearVelocityField(u0=u0, v0=v0);

            maxError = LinearVelocityFieldKnownIssueTests.maxPathError(model,x0,y0,2*86400,360,u0,v0);

            testCase.verifyLessThanOrEqual(maxError, 5e-6);
        end
    end

    methods (Static, Access = private)
        function [x0, y0] = legacyInitialPositions()
            x0 = 4*[-500; -250; 0; 0; 0; 0; 0; 250; 500];
            y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0];
        end

        function maxError = maxPathError(model,x0,y0,T,dt,u0,v0)
            integrator = AdvectionDiffusionIntegrator(model,0);
            [t, x, y] = integrator.particleTrajectories(x0,y0,T,dt);
            [xAnalytical, yAnalytical] = model.particlePath(x0, y0, t, 0, u0, v0);

            maxError = max(hypot(x - xAnalytical, y - yAnalytical), [], "all");
        end
    end

end
