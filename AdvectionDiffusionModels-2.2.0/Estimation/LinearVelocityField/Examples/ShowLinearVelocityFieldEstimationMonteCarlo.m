%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show direct linear velocity field estimation across Monte Carlo trials.
%
% This example generates synthetic trajectories for several strain-only and
% strain-plus-vorticity cases, estimates the direct linear velocity field
% parameters, and plots the recovered parameter distributions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nModels = 3;
referenceLineColor = [0 0 0];
referenceLineWidth = 2;

for iModel = 1:nModels
    if iModel == 1
        % Three different parameters values in a strain only model.
        modelName = 'strain-diffusive';
        model = [ModelParameter.strain];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = 0*zeros(size(sigmaValues));
        thetaValues = [45 0 -45]*pi/180;
        T = 3*86400;
        dt = 4*3600;
    elseif iModel == 2
        % Parameters values for a velocity field with both strain and
        % vorticity, but strain dominated.
        modelName = 'vorticity-strain-diffusive';
        model = [ModelParameter.strain,ModelParameter.vorticity];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -0.5*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    elseif iModel == 3
        % Parameters values for a velocity field with both strain and
        % vorticity, but vortity dominated.
        modelName = 'vorticity-strain-diffusive';
        model = [ModelParameter.strain,ModelParameter.vorticity];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -2*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    end
    
    
    % u0 = 0.1*0;
    % v0 = -0.2*0;
    % sigma = 1e-5;
    % zeta = -3e-6;
    % theta = 30*pi/180;
    % kappa = 0.1;
    % T = round(1/sigma/86400)*86400;
    % T = 2*86400;
    % dt = 3600;
    % model='strain-diffusive';
    % model = [ModelParameter.strain];
    % % model='vorticity-strain-diffusive';
    
    for iParameter = 1:nParameters
        sigma = sigmaValues(iParameter);
        zeta = zetaValues(iParameter);
        kappa = kappaValues(iParameter);
        theta = thetaValues(iParameter);
        plainFigureName = sprintf('%s | case %d of %d', modelName, iParameter, nParameters);
        latexFigureTitle = sprintf('%s, case %d of %d, $\\kappa = %s$, $\\sigma = %s$, $\\theta = %s$, $\\zeta = %s$', modelName, iParameter, nParameters, latexScalar(kappa), latexScalar(sigma), latexAngleDegrees(theta*180/pi), latexScalar(zeta));
        u0 = 0;
        v0 = 0;
        velocityField = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
        integrator = AdvectionDiffusionIntegrator(velocityField,kappa);
        integrator.stepSize = 1800;
        
%         x = linspace(-500,500,5);
%         y = linspace(-500,500,5);
%         [x0,y0] = ndgrid(x,y);
        
        % Particles placed in a "plus" configuration.
        x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
        y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];
        
        
        totalIterations = 100;
        u0Est = zeros(totalIterations,1);
        v0Est = zeros(totalIterations,1);
        sigma_nEst = zeros(totalIterations,1);
        sigma_sEst = zeros(totalIterations,1);
        kappaEst = zeros(totalIterations,1);
        zetaEst = zeros(totalIterations,1);
        
        for i=1:totalIterations
            [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
            
            parameterEstimates = EstimateLinearVelocityFieldParameters( x, y, t, model );
            u0Est(i) = parameterEstimates.u0;
            v0Est(i) = parameterEstimates.v0;
            sigma_nEst(i) = parameterEstimates.sigma_n;
            sigma_sEst(i) = parameterEstimates.sigma_s;
            zetaEst(i) = parameterEstimates.zeta;
            
            [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm] = DecomposeTrajectories(x, y, t, parameterEstimates);
            x_sm = cumtrapz(t,u_sm);
            y_sm = cumtrapz(t,v_sm);
            kappaEst(i) = mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));
        end
        
        
        % Kernel density estimate of the distribution.
        data = cat(2,sigma_nEst,sigma_sEst);
        densityModel = KernelDensityEstimate.fromData(data);
        [density, gridVectors] = densityModel.densityOnGrid();
        
        figure('Name', plainFigureName, 'NumberTitle', 'off', 'Position', [50 50 900 400])
        subplot(2,4,[1 2 5 6])
        contourf(gridVectors{1}, gridVectors{2}, density.');
        hold on
        for i=2:2:10
            rectangle('Position',[-i -i 2*i 2*i]*1e-6, 'Curvature', [1 1]);
        end
        scatter(sigma_nEst, sigma_sEst, 5^2, 0*[1 1 1],'filled');
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 1*[1 1 1],'LineWidth',5);
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 0*[1 1 1],'LineWidth',3);
        axis equal
        xlabel('$\sigma_n$','Interpreter','latex')
        ylabel('$\sigma_s$','Interpreter','latex')
        
        % make the zeros be white
        cmap = colormap;
        cmap(1,:)=1;
        colormap(cmap)
        
        subplot(2,4,3)
        histogram(u0Est,'Normalization','pdf')
        xline(u0, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$u_0$','Interpreter','latex')
        
        subplot(2,4,4)
        histogram(v0Est,'Normalization','pdf')
        xline(v0, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$v_0$','Interpreter','latex')
        
        subplot(2,4,7)
        histogram(zetaEst,'Normalization','pdf')
        xline(zeta, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$\zeta$','Interpreter','latex')
        
        subplot(2,4,8)
        histogram(kappaEst,'Normalization','pdf')
        xline(kappa, 'Color', referenceLineColor, 'LineStyle', '--', 'LineWidth', referenceLineWidth)
        title('$\kappa$','Interpreter','latex')

        sgtitle(latexFigureTitle,'Interpreter','latex')
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
