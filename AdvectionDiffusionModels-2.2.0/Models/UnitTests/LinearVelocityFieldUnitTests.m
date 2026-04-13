classdef LinearVelocityFieldUnitTests < matlab.unittest.TestCase

    methods (Test)
        function derivedParametersMatchSigmaThetaDefinition(testCase)
            sigma = 4e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);

            testCase.verifyEqual(model.sigma_n, sigma*cos(2*theta), "AbsTol", 1e-15);
            testCase.verifyEqual(model.sigma_s, sigma*sin(2*theta), "AbsTol", 1e-15);
        end

        function pointwiseVelocityMatchesAffineFormula(testCase)
            sigma = 4e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);

            x = [-600; 0; 450];
            y = [300; -200; 50];

            uExpected = u0 + 0.5*(model.sigma_n*x + (model.sigma_s - zeta)*y);
            vExpected = v0 + 0.5*((model.sigma_s + zeta)*x - model.sigma_n*y);

            testCase.verifyEqual(model.u(1234,x,y), uExpected, "AbsTol", 1e-12);
            testCase.verifyEqual(model.v(1234,x,y), vExpected, "AbsTol", 1e-12);
        end

        function particlePathMatchesIntegratorForPassingParameterSweep(testCase)
            sigmas = [2e-6 8e-6 2e-5];
            thetas = [pi/4 0 -pi/4];
            zetaFactors = [0 -0.5 -2];
            durations = [3*86400 5*86400 5*86400];
            u0 = 0.08;
            v0 = -0.05;

            [x0, y0] = LinearVelocityFieldUnitTests.legacyInitialPositions();

            for sigma = sigmas
                for theta = thetas
                    for iCase = 1:length(zetaFactors)
                        zeta = zetaFactors(iCase)*sigma;
                        T = durations(iCase);

                        model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
                        maxError = LinearVelocityFieldUnitTests.maxPathError(model,x0,y0,T,360,u0,v0);

                        diagnostic = sprintf('particlePath mismatch for sigma=%g, theta=%g, zeta=%g', sigma, theta, zeta);
                        testCase.verifyLessThanOrEqual(maxError, 5e-6, diagnostic);
                    end
                end
            end
        end

        function particlePathMatchesIntegratorForMatchedFlowWithoutBackgroundVelocity(testCase)
            sigma = 8e-6;
            theta = 0;
            zetas = [sigma -sigma];
            u0 = 0;
            v0 = 0;

            [x0, y0] = LinearVelocityFieldUnitTests.legacyInitialPositions();

            for zeta = zetas
                model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
                maxError = LinearVelocityFieldUnitTests.maxPathError(model,x0,y0,2*86400,360,u0,v0);

                diagnostic = sprintf('Matched-flow particlePath mismatch for sigma=%g, zeta=%g', sigma, zeta);
                testCase.verifyLessThanOrEqual(maxError, 5e-6, diagnostic);
            end
        end

        function particlePathMatchesIntegratorForZeroFieldWithoutBackgroundVelocity(testCase)
            [x0, y0] = LinearVelocityFieldUnitTests.legacyInitialPositions();
            model = LinearVelocityField();
            integrator = AdvectionDiffusionIntegrator(model,0);

            [t, x, y] = integrator.particleTrajectories(x0,y0,2*86400,360);
            [xAnalytical, yAnalytical] = model.particlePath(x0, y0, t, 0, 0, 0);

            xExpected = repmat(x0.', length(t), 1);
            yExpected = repmat(y0.', length(t), 1);

            testCase.verifyEqual(x, xExpected, "AbsTol", 1e-12);
            testCase.verifyEqual(y, yExpected, "AbsTol", 1e-12);
            testCase.verifyEqual(xAnalytical, xExpected, "AbsTol", 1e-12);
            testCase.verifyEqual(yAnalytical, yExpected, "AbsTol", 1e-12);
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
