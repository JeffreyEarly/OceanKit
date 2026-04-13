%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple box
%
% The idea here is that we start a bunch of particles in part of the box
% and after some time, the particles evenly spread to all parts of the box.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 0;
videoPath = 'particles_in_a_box.mp4';

box = SimpleBox();

kappa = 100;
integrator = AdvectionDiffusionIntegrator(box, kappa);

% determine reasonable integration time scales
T = 2 * 86400;
dt = 864;

% place particles on the left-quarter of the domain
x = linspace(min(box.xlim), max(box.xlim) / 4, 20);
y = linspace(min(box.ylim), max(box.ylim), 10);
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

for iTime = iTime0:length(t)
    clf
    
    box.plotBounds();
    hold on
    axis equal
    xticks([])
    yticks([])
    axis off
    
    scatter(x(iTime,:) / box.visualScale, y(iTime,:) / box.visualScale, 8^2, 'r', 'filled')
    
    if showEndStateOnly == 1
        continue;
    end

    drawnow
    writeVideo(videoWriter, getframe(figureHandle));
end

if ~showEndStateOnly
    close(videoWriter);
end
