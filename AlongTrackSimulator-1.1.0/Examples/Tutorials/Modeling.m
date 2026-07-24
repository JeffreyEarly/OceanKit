%% Tutorial Metadata
% Title: Modeling with satellite observations
% Slug: modeling-with-satellite-observations
% Description: Model a mesoscale eddy with WaveVortexModel and compare how current altimetry missions sample it.
% NavOrder: 1

%% Initialize the model domain
% We start with a barotropic quasi-geostrophic domain large enough to watch
% one mesoscale eddy translate across the basin while the active missions
% sample it along their native ground tracks.

Lx = 2000e3;
Ly = 1000e3;
Nx = 2*256;
Ny = 2*128;
latitude = 24;

wvt = WVTransformBarotropicQG([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

%% Add a Gaussian eddy
% The initial sea-surface height anomaly is
%
% $$
% \eta(x,y)=A\exp\left(-\frac{(x-x_0)^2 + (y-y_0)^2}{L^2}\right),
% $$
%
% so the tutorial begins from one isolated anticyclonic feature.

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x, y) A*exp(-((x-x0).^2 + (y-y0).^2)/L^2), shouldRemoveMeanPressure=true);

figure(Position=[100 100 760 320])
pcolor(wvt.x/1e3, wvt.y/1e3, 100*(wvt.ssh).'), shading interp
axis equal tight
xlabel("x (km)")
ylabel("y (km)")
title("Initial sea-surface height anomaly")
colorbar
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("initial-eddy-state", Caption="The Gaussian eddy starts as one compact sea-surface height anomaly centered away from the boundaries."); end

%% Add the background dynamics
% We add adaptive damping and beta-plane advection before launching the
% model integration.

wvt.addForcing(WVAdaptiveDamping(wvt));
wvt.addForcing(WVBetaPlanePVAdvection(wvt));
model = WVModel(wvt);

%% Add the active altimetry missions
% The along-track output group samples the model with the geometry of each
% currently operating mission.

netcdfPath = string(tempname) + ".nc";
cleanupNetCDF = onCleanup(@() deleteTemporaryFile(netcdfPath)); %#ok<NASGU>
outputFile = model.createNetCDFFileForModelOutput(netcdfPath, outputInterval=86400, shouldOverwriteExisting=true);
model.addNetCDFOutputVariables("ssh", "zeta_z");

ats = AlongTrackSimulator();
currentMissions = ats.currentMissions;
ats.addMissionsToOutputFile(outputFile, missionNames=currentMissions);

%% Integrate for one year
% This gives the eddy time to translate westward while the repeat and
% drifting missions build up different spatial sampling patterns.

model.integrateToTime(365*86400);
outputFile.closeNetCDFFile();

figure(Position=[100 100 760 320])
pcolor(wvt.x/1e3, wvt.y/1e3, 100*(wvt.ssh).'), shading interp
axis equal tight
xlabel("x (km)")
ylabel("y (km)")
title("Sea-surface height after one year")
colorbar
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("one-year-eddy-state", Caption="After one year the eddy has propagated westward and broadened under the model dynamics."); end

%% Render the satellite sampling movie
% The movie combines the full field with three observation views: the
% reference mission, the reference-plus-interleaved pair, and all active
% missions together.

movieSettings = modelingMovieSettings(currentMissions);
if exist("tutorialMovieCapture", "var") && isa(tutorialMovieCapture, "function_handle"), tutorialMovieCapture( ...
        "QGMonopoleSatelliteObservations.mp4", ...
        Caption="Daily snapshots of the full field are shown with the reference mission, the reference-plus-interleaved pair, and all active missions using a four-day fading window around each frame.", ...
        Poster="t-001.png", ...
        Settings=movieSettings, ...
        Build=@(moviePath, posterPath) buildModelingMovie(netcdfPath, moviePath, posterPath, movieSettings)); end

% The upper-left panel shows the full model field.
%
% The upper-right panel shows the reference mission, Sentinel-6A, with a
% four-day fade before and after each frame so the repeat-track pattern is
% visible.
%
% The lower-left panel adds the interleaved Jason-3 orbit to show how the
% staggered pair fills gaps between successive Sentinel-6A passes.
%
% The lower-right panel overlays all active missions to show how the full
% constellation samples the eddy.

function settings = modelingMovieSettings(allMissions)
settings = struct( ...
    "AllMissions", reshape(string(allMissions), 1, []), ...
    "ReferenceMissions", "s6a", ...
    "ReferencePlusInterleavedMissions", ["s6a", "j3n"], ...
    "WindowDays", 4, ...
    "FrameRate", 12, ...
    "Quality", 60, ...
    "ScatterSize", 3^2, ...
    "ColorLimits", [-3 15], ...
    "FigurePosition", [50 50 860 420], ...
    "Resolution", 150);
end

function buildModelingMovie(netcdfPath, moviePath, posterPath, settings)
[wvt, ncfile] = WVTransform.waveVortexTransformFromFile(netcdfPath);
cleanupNetCDF = onCleanup(@() closeNetCDFIfPossible(ncfile)); %#ok<NASGU>
missionTracks = readMissionTracks(ncfile, settings.AllMissions);
tGridded = ncfile.readVariables("wave-vortex/t");
movieTimeIndices = find(tGridded <= (tGridded(end) - settings.WindowDays*86400));
if isempty(movieTimeIndices)
    movieTimeIndices = 1:numel(tGridded);
end

fig = figure(Units="points", Position=settings.FigurePosition, Visible="off", Color="w");
set(fig, "PaperPositionMode", "auto")
cleanupFigure = onCleanup(@() closeIfValid(fig)); %#ok<NASGU>

videoWriter = VideoWriter(moviePath, "MPEG-4");
videoWriter.FrameRate = settings.FrameRate;
videoWriter.Quality = settings.Quality;
open(videoWriter);
cleanupVideoWriter = onCleanup(@() closeVideoWriterIfOpen(videoWriter)); %#ok<NASGU>

frameImagePath = string(tempname) + ".png";
cleanupFrameImage = onCleanup(@() deleteTemporaryFile(frameImagePath)); %#ok<NASGU>
didWritePoster = false;

for iTime = reshape(movieTimeIndices, 1, [])
    wvt.initFromNetCDFFile(ncfile, iTime=iTime)
    renderModelingMovieFrame(fig, wvt, missionTracks, settings);
    exportgraphics(fig, frameImagePath, Resolution=settings.Resolution);
    if ~didWritePoster && posterPath ~= ""
        copyfile(frameImagePath, posterPath, "f");
        didWritePoster = true;
    end
    writeVideo(videoWriter, imread(frameImagePath));
end

close(videoWriter)
clear cleanupVideoWriter
closeIfValid(fig)
clear cleanupFigure
closeNetCDFIfPossible(ncfile)
clear cleanupNetCDF
end

function missionTracks = readMissionTracks(ncfile, missionNames)
missionNames = reshape(string(missionNames), [], 1);
missionTracks = repmat(struct("MissionName", "", "t", [], "x", [], "y", [], "ssh", []), numel(missionNames), 1);
for iMission = 1:numel(missionNames)
    group = ncfile.groupWithName(missionNames(iMission));
    [missionTracks(iMission).t, missionTracks(iMission).x, missionTracks(iMission).y, missionTracks(iMission).ssh] = ...
        group.readVariables("t", "track_x", "track_y", "ssh");
    missionTracks(iMission).MissionName = missionNames(iMission);
end
end

function renderModelingMovieFrame(fig, wvt, missionTracks, settings)
clf(fig)
tileLayout = tiledlayout(fig, 2, 2, TileSpacing="tight", Padding="compact");
title(tileLayout, sprintf("day %.2f", wvt.t/86400))
colorMap = parula(256);

ax = nexttile(tileLayout, 1);
renderFullFieldPanel(ax, wvt, colorMap, settings);

ax = nexttile(tileLayout, 2);
renderObservationPanel(ax, wvt, missionTracks, settings.ReferenceMissions, colorMap, settings, "reference mission", false, false);

ax = nexttile(tileLayout, 3);
renderObservationPanel(ax, wvt, missionTracks, settings.ReferencePlusInterleavedMissions, colorMap, settings, "reference + interleaved", true, true);

ax = nexttile(tileLayout, 4);
renderObservationPanel(ax, wvt, missionTracks, settings.AllMissions, colorMap, settings, "all active missions", true, false);
end

function renderFullFieldPanel(ax, wvt, colorMap, settings)
pcolor(ax, wvt.x/1e3, wvt.y/1e3, 100*(wvt.ssh).')
shading(ax, "interp")
colormap(ax, colorMap)
clim(ax, settings.ColorLimits)
xlim(ax, [min(wvt.x) max(wvt.x)]/1e3)
ylim(ax, [min(wvt.y) max(wvt.y)]/1e3)
axis(ax, "equal", "tight")
ax.XTick = [];
ylabel(ax, "km")
title(ax, "full field")
end

function renderObservationPanel(ax, wvt, missionTracks, missionNames, colorMap, settings, panelTitle, showXAxis, showYAxis)
missionNames = reshape(string(missionNames), 1, []);
windowSeconds = settings.WindowDays*86400;
hold(ax, "on")

for iMission = 1:numel(missionTracks)
    if ~ismember(missionTracks(iMission).MissionName, missionNames)
        continue;
    end

    indicesInRange = abs(missionTracks(iMission).t - wvt.t) < windowSeconds;
    if ~any(indicesInRange)
        continue;
    end

    x = missionTracks(iMission).x(indicesInRange);
    y = missionTracks(iMission).y(indicesInRange);
    t = missionTracks(iMission).t(indicesInRange);
    ssh = missionTracks(iMission).ssh(indicesInRange);
    alpha = max(0, 1 - abs(t - wvt.t)/windowSeconds);
    markerHandle = scatter(ax, x/1e3, y/1e3, settings.ScatterSize, 100*ssh, "filled");
    markerHandle.AlphaDataMapping = "none";
    markerHandle.MarkerFaceAlpha = "flat";
    markerHandle.MarkerEdgeAlpha = "flat";
    markerHandle.AlphaData = alpha;
end

colormap(ax, colorMap)
clim(ax, settings.ColorLimits)
xlim(ax, [min(wvt.x) max(wvt.x)]/1e3)
ylim(ax, [min(wvt.y) max(wvt.y)]/1e3)
axis(ax, "equal", "tight")
if showXAxis
    xlabel(ax, "km")
else
    ax.XTick = [];
end
if showYAxis
    ylabel(ax, "km")
else
    ax.YTick = [];
end
title(ax, panelTitle)
end

function closeIfValid(fig)
if isgraphics(fig)
    close(fig)
end
end

function closeVideoWriterIfOpen(videoWriter)
try
    close(videoWriter)
catch
end
end

function closeNetCDFIfPossible(ncfile)
if ismethod(ncfile, "close")
    try
        ncfile.close();
    catch
    end
end
end

function deleteTemporaryFile(filePath)
if isfile(filePath)
    delete(filePath)
end
end
