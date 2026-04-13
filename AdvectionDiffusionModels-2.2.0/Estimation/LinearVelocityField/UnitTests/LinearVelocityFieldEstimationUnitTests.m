classdef LinearVelocityFieldEstimationUnitTests < matlab.unittest.TestCase

    methods (Test)
        function estimateConstantLinearVelocityFieldParameters(testCase)
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model,0);

            [x0, y0] = ndgrid([-600 0 600],[-400 0 400]);
            [t, x, y] = integrator.particleTrajectories(x0(:),y0(:),12*3600,900);

            parameters = EstimateLinearVelocityFieldParameters(x, y, t, [ModelParameter.u0v0, ModelParameter.strain, ModelParameter.vorticity], 1);

            testCase.verifyEqual(parameters.u0, u0, "AbsTol", 5e-6);
            testCase.verifyEqual(parameters.v0, v0, "AbsTol", 5e-6);
            testCase.verifyEqual(parameters.sigma_n, model.sigma_n, "AbsTol", 5e-10);
            testCase.verifyEqual(parameters.sigma_s, model.sigma_s, "AbsTol", 5e-10);
            testCase.verifyEqual(parameters.zeta, zeta, "AbsTol", 5e-10);
        end

        function estimateTimeVaryingBackgroundVelocityWithSplines(testCase)
            t = linspace(0,12*3600,49).';
            T = t(end);

            u0 = 0.05 + 0.03*(t/T);
            v0 = -0.04 + 0.02*(t/T);

            xCenter = 0.05*t + 0.03*(t.^2)/(2*T);
            yCenter = -0.04*t + 0.02*(t.^2)/(2*T);

            [xOffset, yOffset] = ndgrid([-300 0 300],[-200 0 200]);
            x = xCenter + reshape(xOffset,1,[]);
            y = yCenter + reshape(yOffset,1,[]);

            parameters = EstimateLinearVelocityFieldParameters(x, y, t, ModelParameter.u0v0, 4);

            testCase.verifySize(parameters.u0,size(t));
            testCase.verifySize(parameters.v0,size(t));
            testCase.verifyLessThanOrEqual(max(abs(parameters.u0 - u0)), 5e-7);
            testCase.verifyLessThanOrEqual(max(abs(parameters.v0 - v0)), 5e-7);
            testCase.verifyEqual(parameters.sigma_n, zeros(size(t)), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.sigma_s, zeros(size(t)), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.zeta, zeros(size(t)), "AbsTol", 1e-12);
        end

        function estimateSolutionLikelihoodFromBootstrapsMatchesDirectGaussianReference(testCase)
            bootstraps = struct( ...
                "u0", [0.10 0.12 0.15 0.13; 0.20 0.18 0.22 0.25], ...
                "v0", [-0.06 -0.04 -0.05 -0.03; -0.02 -0.01 -0.03 0.00], ...
                "u1", [0.02 0.03 0.01 0.04; -0.01 0.00 0.02 0.01], ...
                "v1", [-0.03 -0.01 -0.02 0.00; 0.01 0.00 0.02 0.03], ...
                "sigma_n", [1.0 1.3 1.1 1.4; 0.8 0.9 1.0 1.2] * 1e-6, ...
                "sigma_s", [-0.5 -0.2 -0.1 0.1; 0.2 0.3 0.4 0.5] * 1e-6, ...
                "zeta", [-2.0 -1.7 -1.9 -1.5; -0.7 -0.6 -0.4 -0.3] * 1e-6, ...
                "delta", [0.6 0.8 0.7 0.9; 0.1 0.3 0.2 0.4] * 1e-6);
            parameters = [ModelParameter.u0v0, ModelParameter.u1v1, ModelParameter.strain, ModelParameter.vorticity, ModelParameter.divergence];

            actual = EstimateSolutionLikelihoodFromBootstraps(bootstraps, parameters, 1);
            expected = LinearVelocityFieldEstimationUnitTests.directLikelihoodFromBootstraps(bootstraps, parameters, 1);

            testCase.verifyEqual(actual, expected, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(actual)))
        end
    end

    methods (Static, Access = private)
        function jointLikelihood = directLikelihoodFromBootstraps(bootstraps, parameters, stride)
            shouldEstimateU0V0 = any(parameters == ModelParameter.u0v0);
            shouldEstimateU1V1 = any(parameters == ModelParameter.u1v1);
            shouldEstimateStrain = any(parameters == ModelParameter.strain);
            shouldEstimateVorticity = any(parameters == ModelParameter.vorticity);
            shouldEstimateDivergence = any(parameters == ModelParameter.divergence);

            nT = size(bootstraps.u0, 1);
            nBootstraps = size(bootstraps.u0, 2);
            jointLikelihood = zeros(length(1:stride:nT), nBootstraps);

            iOutput = 1;
            for iTime = 1:stride:nT
                if shouldEstimateU0V0
                    jointLikelihood(iOutput, :) = jointLikelihood(iOutput, :) + LinearVelocityFieldEstimationUnitTests.directLogDensity2D( ...
                        bootstraps.u0(iTime, :).', bootstraps.v0(iTime, :).');
                end

                if shouldEstimateU1V1
                    jointLikelihood(iOutput, :) = jointLikelihood(iOutput, :) + LinearVelocityFieldEstimationUnitTests.directLogDensity2D( ...
                        bootstraps.u1(iTime, :).', bootstraps.v1(iTime, :).');
                end

                if shouldEstimateStrain
                    jointLikelihood(iOutput, :) = jointLikelihood(iOutput, :) + LinearVelocityFieldEstimationUnitTests.directLogDensity2D( ...
                        bootstraps.sigma_n(iTime, :).', bootstraps.sigma_s(iTime, :).');
                end

                if shouldEstimateVorticity
                    jointLikelihood(iOutput, :) = jointLikelihood(iOutput, :) + LinearVelocityFieldEstimationUnitTests.directLogDensity1D( ...
                        bootstraps.zeta(iTime, :).');
                end

                if shouldEstimateDivergence
                    jointLikelihood(iOutput, :) = jointLikelihood(iOutput, :) + LinearVelocityFieldEstimationUnitTests.directLogDensity1D( ...
                        bootstraps.delta(iTime, :).');
                end

                iOutput = iOutput + 1;
            end

            jointLikelihood = sum(jointLikelihood, 1);
        end

        function score = directLogDensity1D(values)
            model = KernelDensityEstimate.fromData(values);
            pointDensity = LinearVelocityFieldEstimationUnitTests.directGaussianDensity1D(model.data, model.bandwidth, values);
            score = log10(LinearVelocityFieldEstimationUnitTests.safeDensity(pointDensity)).';
        end

        function score = directLogDensity2D(xValues, yValues)
            model = KernelDensityEstimate.fromData([xValues, yValues]);
            pointDensity = LinearVelocityFieldEstimationUnitTests.directGaussianDensity2D(model.data, model.bandwidth, [xValues, yValues]);
            score = log10(LinearVelocityFieldEstimationUnitTests.safeDensity(pointDensity)).';
        end

        function density = safeDensity(values)
            density = reshape(values, [], 1);
            density(~isfinite(density) | density <= 0) = realmin("double");
        end

        function density = directGaussianDensity1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            density = mean(exp(-0.5 * ((queryPoints - data.')/bandwidth).^2), 2)/(bandwidth * sqrt(2*pi));
        end

        function density = directGaussianDensity2D(data, bandwidth, queryPoints)
            density = mean(exp(-0.5 * (((queryPoints(:, 1) - data(:, 1).')/bandwidth(1)).^2 + ((queryPoints(:, 2) - data(:, 2).')/bandwidth(2)).^2)), 2) ...
                /(2*pi * bandwidth(1) * bandwidth(2));
        end
    end

end
