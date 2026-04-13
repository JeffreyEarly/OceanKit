%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Divergence/Convergence box
%
% Here we start with a bunch of particles evenly spread throughout the
% domain, but in this case, there's a weak convergence in half the box and
% a weak divergence in the other half of the box.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 1;
videoPath = 'divergence_convergence.mp4';

model = DivergenceBox();

kappa = 100;
integrator = AdvectionDiffusionIntegrator(model, kappa);

% determine reasonable integration time scales
T = 20*300;
dt = 20;

% place particles on the left-quarter of the domain
x = linspace(min(model.xlim),max(model.xlim),20);
y = linspace(min(model.ylim),max(model.ylim),10);
[x0,y0] = ndgrid(x,y);

n = floor(T/dt/3);
model.U = 0.33;
[t1, x1, y1] = integrator.particleTrajectories(x0, y0, n * dt, dt);
model.U = 0.66;
[t2, x2, y2] = integrator.particleTrajectories(x1(end,:), y1(end,:), n * dt, dt);
model.U = 1.0;
[t3, x3, y3] = integrator.particleTrajectories(x2(end,:), y2(end,:), n * dt, dt);

t = cat(1,t1,t2,t3);
x = cat(1,x1,x2,x3);
y = cat(1,y1,y2,y3);

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
    
    model.plotBounds();
    hold on
    model.plotVelocityField
    axis equal
    xticks([])
    yticks([])
    title([])
    axis off
    
    scatter(x(iTime,:) / model.visualScale, y(iTime,:) / model.visualScale, 8^2, 'r', 'filled')
    
    if iTime < n
        text(0, model.Ly / model.visualScale - 0.5, 'Weak divergence')
    elseif iTime < 2 * n
        text(0, model.Ly / model.visualScale - 0.5, 'Moderate divergence')
    else
        text(0, model.Ly / model.visualScale - 0.5, 'Strong divergence')
    end
    
    if showEndStateOnly == 1
        continue;
    end

    drawnow
    writeVideo(videoWriter, getframe(figureHandle));
end

if ~showEndStateOnly
    close(videoWriter);
end
