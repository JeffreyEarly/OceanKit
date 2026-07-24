[wvt,ncfile] = WVTransform.waveVortexTransformFromFile('QGMonopoleWithAlongTrack.nc');
t_gridded = ncfile.readVariables("wave-vortex/t");

%%
allMissions = ["alg","c2n","h2b","j3n","s3a","s3b","s6a","swon"];
tracks = cell(length(allMissions),1);
for iMission=1:length(allMissions)
    group = ncfile.groupWithName(allMissions(iMission));
    [tracks{iMission}.t, tracks{iMission}.x, tracks{iMission}.y, tracks{iMission}.ssh] = group.readVariables("t","track_x","track_y","ssh");
end

%%

panel1Missions = "s6a";
panel2Missions = ["s6a","j3n"];
panel3Missions = allMissions;

T = 4*86400;
scatSz = 3^2;
colorlimits = [-3 15];
t_windowed = t_gridded;
t_windowed(t_windowed > (t_windowed(end)-T)) = [];

fig1 = figure('Units', 'points', 'Position', [50 50 860 420],Visible='off');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

tl = tiledlayout(2,2,TileSpacing="tight");
for iTime = 1:length(t_windowed)
    wvt.initFromNetCDFFile(ncfile,iTime=iTime)
    title(tl,sprintf("day %.2f",wvt.t/86400))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % full field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile(tl,1), cla
    % griddedAxes = axes;
    pcolor(wvt.x/1e3,wvt.y/1e3,100*(wvt.ssh).');
    ax = gca;
    shading(ax,"interp")
    cm = colormap(ax,"parula");
    % colormap((1-alpha)*cm + alpha*ones(size(cm,1),3));
    clim(ax,colorlimits);
    ax.XTick = [];
    ax.XTickLabel = [];
    ylabel('km')
    xlim(ax,[min(wvt.x) max(wvt.x)/1e3])
    ylim(ax,[min(wvt.y) max(wvt.y)/1e3])
    axis(ax,'equal','tight');
    title('full-field')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reference mission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile(tl,2), cla
    [~,missionIndices] = ismember(panel1Missions,allMissions);
    for iIndex=1:length(missionIndices)
        iMission = missionIndices(iIndex);
        indicesInRange = abs(tracks{iMission}.t-wvt.t) < T;
        x = tracks{iMission}.x(indicesInRange);
        y = tracks{iMission}.y(indicesInRange);
        t = tracks{iMission}.t(indicesInRange);
        ssh = tracks{iMission}.ssh(indicesInRange);
        alpha = max(0, 1 - abs(t - wvt.t)/T);
        h = scatter(x/1e3, y/1e3, scatSz, 100*ssh, 'filled' ); hold on;

        h.AlphaDataMapping = 'none';
        h.MarkerFaceAlpha  = 'flat';    % per‐marker face alpha
        h.MarkerEdgeAlpha  = 'flat';    % per‐marker edge alpha
        h.AlphaData        = alpha;
    end
    colormap(cm);
    clim(colorlimits);

    xlim([min(wvt.x) max(wvt.x)/1e3])
    ylim([min(wvt.y) max(wvt.y)/1e3])
    axis('equal','tight');
    ax = gca;
    ax.XTick = []; ax.YTick = [];
    title("reference mission")

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reference + interleaved mission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile(tl,3), cla
    [~,missionIndices] = ismember(panel2Missions,allMissions);
    for iIndex=1:length(missionIndices)
        iMission = missionIndices(iIndex);
        indicesInRange = abs(tracks{iMission}.t-wvt.t) < T;
        x = tracks{iMission}.x(indicesInRange);
        y = tracks{iMission}.y(indicesInRange);
        t = tracks{iMission}.t(indicesInRange);
        ssh = tracks{iMission}.ssh(indicesInRange);
        alpha = max(0, 1 - abs(t - wvt.t)/T);
        h = scatter(x/1e3, y/1e3, scatSz, 100*ssh, 'filled' ); hold on;

        h.AlphaDataMapping = 'none';
        h.MarkerFaceAlpha  = 'flat';    % per‐marker face alpha
        h.MarkerEdgeAlpha  = 'flat';    % per‐marker edge alpha
        h.AlphaData        = alpha;
    end
    colormap(cm);
    clim(colorlimits);

    xlim([min(wvt.x) max(wvt.x)/1e3])
    ylim([min(wvt.y) max(wvt.y)/1e3])
    axis('equal','tight');
    % ax.XTick = []; ax.YTick = [];
    xlabel('km')
    ylabel('km')
    title("reference + interleaved mission")


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % all missions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile(tl,4), cla
    [~,missionIndices] = ismember(panel3Missions,allMissions);
    for iIndex=1:length(missionIndices)
        iMission = missionIndices(iIndex);
        indicesInRange = abs(tracks{iMission}.t-wvt.t) < T;
        x = tracks{iMission}.x(indicesInRange);
        y = tracks{iMission}.y(indicesInRange);
        t = tracks{iMission}.t(indicesInRange);
        ssh = tracks{iMission}.ssh(indicesInRange);
        alpha = max(0, 1 - abs(t - wvt.t)/T);
        h = scatter(x/1e3, y/1e3, scatSz, 100*ssh, 'filled' ); hold on;

        h.AlphaDataMapping = 'none';
        h.MarkerFaceAlpha  = 'flat';    % per‐marker face alpha
        h.MarkerEdgeAlpha  = 'flat';    % per‐marker edge alpha
        h.AlphaData        = alpha;
    end
    colormap(cm);
    clim(colorlimits);

    xlim([min(wvt.x) max(wvt.x)/1e3])
    ylim([min(wvt.y) max(wvt.y)/1e3])
    axis('equal','tight');
    ax=gca;
    ax.YTick = [];
    xlabel('km')
    title("all missions")

    exportgraphics(fig1,sprintf('./movie-figures/t-%03d.png',iTime),Resolution=300)
end