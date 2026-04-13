%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show gridded-streamfunction bootstrap estimation across synthetic cases.
%
% This example mirrors the synthetic case families and figure layout from
% `ShowLinearVelocityFieldEstimationMonteCarlo`, but replaces the direct
% linear estimator with `GriddedStreamfunctionBootstrap`.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nModels = 3;
nBootstraps = 100;
referenceLineColor = [0 0 0];
referenceLineWidth = 2;

for iModel = 1:nModels
    if iModel == 1
        modelName = 'strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = 0*zeros(size(sigmaValues));
        thetaValues = [45 0 -45]*pi/180;
        T = 3*86400;
        dt = 4*3600;
    elseif iModel == 2
        modelName = 'vorticity-strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -0.5*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    else
        modelName = 'vorticity-strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -2*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    end

    for iParameter = 1:nParameters
        sigma = sigmaValues(iParameter);
        zeta = zetaValues(iParameter);
        kappa = kappaValues(iParameter);
        theta = thetaValues(iParameter);
        plainFigureName = sprintf('%s | case %d of %d', modelName, iParameter, nParameters);
        latexFigureTitle = sprintf([ ...
            '%s, case %d of %d, $\\kappa = %s$, $\\sigma = %s$, ' ...
            '$\\theta = %s$, $\\zeta = %s$'], ...
            modelName, iParameter, nParameters, latexScalar(kappa), ...
            latexScalar(sigma), latexAngleDegrees(theta*180/pi), latexScalar(zeta));

        u0 = 0;
        v0 = 0;
        velocityField = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
        integrator = AdvectionDiffusionIntegrator(velocityField, kappa);
        integrator.stepSize = 1800;

        x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
        y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0];

        [t, x, y] = integrator.particleTrajectories(x0, y0, T, dt);
        trajectoryCell = cell(size(x, 2), 1);
        for iDrifter = 1:size(x, 2)
            trajectoryCell{iDrifter} = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
        end
        trajectories = vertcat(trajectoryCell{:});

        bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
            trajectories, nBootstraps=nBootstraps, randomSeed=100*iModel + iParameter);

        sigmaNMean = mean(bootstrap.summary.sigma_n, 1);
        sigmaSMean = mean(bootstrap.summary.sigma_s, 1);
        uCenterMean = mean(bootstrap.summary.uCenter, 1);
        vCenterMean = mean(bootstrap.summary.vCenter, 1);
        zetaMean = mean(bootstrap.summary.zeta, 1);
        kappaEstimate = bootstrap.bootstrapKappa;

        data = [sigmaNMean(:), sigmaSMean(:)];
        densityModel = KernelDensityEstimate.fromData(data);
        [density, gridVectors] = densityModel.densityOnGrid();

        figure('Name', plainFigureName, 'NumberTitle', 'off', 'Position', [50 50 900 400])
        subplot(2,4,[1 2 5 6])
        contourf(gridVectors{1}, gridVectors{2}, density.');
        hold on
        for iRadius = 2:2:10
            rectangle('Position', [-iRadius -iRadius 2*iRadius 2*iRadius]*1e-6, 'Curvature', [1 1]);
        end
        scatter(sigmaNMean, sigmaSMean, 5^2, 0*[1 1 1], 'filled');
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 1*[1 1 1], 'LineWidth', 5);
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 0*[1 1 1], 'LineWidth', 3);
        axis equal
        xlabel('$\sigma_n$', 'Interpreter', 'latex')
        ylabel('$\sigma_s$', 'Interpreter', 'latex')
        cmap = colormap;
        cmap(1,:) = 1;
        colormap(cmap)

        subplot(2,4,3)
        histogram(uCenterMean, 'Normalization', 'pdf')
        xline(u0, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$u_c$', 'Interpreter', 'latex')

        subplot(2,4,4)
        histogram(vCenterMean, 'Normalization', 'pdf')
        xline(v0, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$v_c$', 'Interpreter', 'latex')

        subplot(2,4,7)
        histogram(zetaMean, 'Normalization', 'pdf')
        xline(zeta, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$\zeta$', 'Interpreter', 'latex')

        subplot(2,4,8)
        histogram(kappaEstimate, 'Normalization', 'pdf')
        xline(kappa, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$\kappa$', 'Interpreter', 'latex')

        sgtitle(latexFigureTitle, 'Interpreter', 'latex')
    end
end

function valueString = latexScalar(value)
if value == 0
    valueString = '0';
    return
end

if abs(value) >= 1e-2 && abs(value) < 1e2
    valueString = sprintf('%.2g', value);
    return
end

exponent = floor(log10(abs(value)));
mantissa = value/10^exponent;
if abs(mantissa - round(mantissa)) < 1e-12
    mantissaString = sprintf('%d', round(mantissa));
else
    mantissaString = sprintf('%.2g', mantissa);
end

valueString = sprintf('%s \\times 10^{%d}', mantissaString, exponent);
end

function angleString = latexAngleDegrees(angleDegrees)
if abs(angleDegrees - round(angleDegrees)) < 1e-12
    angleString = sprintf('%d^\\circ', round(angleDegrees));
else
    angleString = sprintf('%.1f^\\circ', angleDegrees);
end
end
