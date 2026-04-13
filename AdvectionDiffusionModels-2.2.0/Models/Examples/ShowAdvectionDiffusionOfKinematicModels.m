shouldSaveImages = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Meandering Jet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jet = MeanderingJet();
kappa = 1e3;
integrator = AdvectionDiffusionIntegrator(jet, kappa);

% determine reasonable integration time scales
T = 5 * jet.Lx / jet.U;
dt = 864;

% place particles throughout the valid domain
x = linspace(min(jet.xlim), max(jet.xlim), 6);
y = linspace(min(jet.ylim), max(jet.ylim), 6);
[x0, y0] = ndgrid(x, y);

[t, x, y] = integrator.particleTrajectories(x0, y0, T, dt);

figureHandle = figure;
jet.plotVelocityField();
hold on
jet.plotTrajectories(x, y)
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/trajectories_jet.png')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Translating Eddy (Gaussian)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eddy = TranslatingGaussian();
% eddy.cx = 0;

kappa = 1e2;
integrator = AdvectionDiffusionIntegrator(eddy, kappa);

% place particles throughout the valid domain
x0 = zeros(1,10);
y0 = linspace(min(eddy.yVisualLimits), max(eddy.yVisualLimits), 10);

T = 200*86400;
dt = 864;

[t, x, y] = integrator.particleTrajectories(x0, y0, T, dt);

figure
eddy.plotTrajectories(x, y)

% cylinder = CylinderFlow();
