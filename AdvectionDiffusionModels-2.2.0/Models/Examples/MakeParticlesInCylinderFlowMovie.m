%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cylinder Flow
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 1;
videoPath = 'particles_in_cylinder_flow.mp4';

model = CylinderFlow();

% let's make the top and bottom non-periodic.
model.xVisualLimits = 2 * model.xVisualLimits;
model.xlim = model.xVisualLimits;
model.ylim = model.yVisualLimits;
model.isYPeriodic = false;

kappa = 1e3;
integrator = AdvectionDiffusionIntegrator(model, kappa);

% determine reasonable integration time scales
T = 4 * 86400;
dt = 864;

sigma_dx = sqrt(2 * kappa * dt);
fprintf('Typical particle step is %.2f km\n', sigma_dx * 1e-3)

% place particles on the left-quarter of the domain
x = linspace(min(model.xlim), max(model.xlim), 20);
y = linspace(min(model.ylim), max(model.ylim), 10);
[x0, y0] = ndgrid(x, y);

[t, x, y] = integrator.particleTrajectories(x0, y0, T, dt);

figureHandle = figure('Position', [50 50 680 400]);
set(figureHandle, 'PaperPositionMode', 'auto')
set(figureHandle, 'Color', 'w');


if showEndStateOnly == 1
    iTime0 = length(t);
else
    iTime0 = 1;
    videoWriter = VideoWriter(videoPath, 'MPEG-4');
    videoWriter.FrameRate = 30;
    open(videoWriter);
end

tailLength = 50;
for iTime = iTime0:length(t)
    clf
    
    model.plotBounds('LineWidth', 8, 'Color', 'black');
    hold on
    axis equal
    xticks([])
    yticks([])
    axis off
    
    startIndex = iTime - tailLength;
    if startIndex < 1
        startIndex = 1;
    end
    range = startIndex:iTime;
    
    model.plotVelocityField();
    hold on
    model.plotTrajectories(x(range,:), y(range,:), 'LineWidth', 1.5)
    
    if showEndStateOnly == 1
        continue;
    end

    drawnow
    writeVideo(videoWriter, getframe(figureHandle));
end

if ~showEndStateOnly
    close(videoWriter);
end
