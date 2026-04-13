classdef SecondMomentMethodUnitTests < matlab.unittest.TestCase

    methods (Test)
        function fitSecondMomentToLinearizedEllipseModelRecoversAnalyticalMoments(testCase)
            sigma = 2e-6;
            theta = pi/6;
            kappa = 1;
            Mxx0 = 1.1e5;
            Myy0 = 0.9e5;
            Mxy0 = 0.2e5;
            t = linspace(0,10*86400,30).';

            [Mxx, Myy, Mxy] = MomentTensorEvolutionInStrainVorticityField(Mxx0, Myy0, Mxy0, t, 0, sigma, theta, kappa);
            [parameters, error] = FitSecondMomentToLinearizedEllipseModel(Mxx, Myy, Mxy, t, 'strain-diffusive');

            testCase.verifyEqual(error, 0, "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.kappa, kappa, "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.sigma*cos(2*parameters.theta), sigma*cos(2*theta), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.sigma*sin(2*parameters.theta), sigma*sin(2*theta), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.zeta, 0, "AbsTol", 1e-12);
        end

        function fitTrajectoriesToLinearizedEllipseModelRecoversZeroDiffusionTrajectories(testCase)
            sigma = 4e-6;
            theta = pi/9;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=0);
            integrator = AdvectionDiffusionIntegrator(model,0);
            [x0, y0] = SecondMomentMethodUnitTests.legacyInitialPositions();
            [t, x, y] = integrator.particleTrajectories(x0,y0,3*86400,1800);

            [parameters, error] = FitTrajectoriesToLinearizedEllipseModel(x, y, t, 'strain-diffusive');

            testCase.verifyEqual(error, 0, "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.kappa, 0, "AbsTol", 5e-12);
            testCase.verifyEqual(parameters.sigma*cos(2*parameters.theta), model.sigma_n, "AbsTol", 5e-12);
            testCase.verifyEqual(parameters.sigma*sin(2*parameters.theta), model.sigma_s, "AbsTol", 5e-12);
            testCase.verifyEqual(parameters.zeta, 0, "AbsTol", 1e-12);
        end

        function fitTrajectoriesToLinearVelocityFieldAcceptsModernSplineFactories(testCase)
            sigma = 4e-6;
            theta = pi/9;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=0);
            integrator = AdvectionDiffusionIntegrator(model,0);
            [x0, y0] = SecondMomentMethodUnitTests.legacyInitialPositions();
            [t, x, y] = integrator.particleTrajectories(x0,y0,3*86400,1800);

            testCase.verifyWarningFree(@() FitTrajectoriesToLinearVelocityField(x, y, t, 'strain-diffusive', 6, 3))
        end
    end

    methods (Static, Access = private)
        function [x0, y0] = legacyInitialPositions()
            x0 = 4*[-500; -250; 0; 0; 0; 0; 0; 250; 500];
            y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0];
        end
    end

end
