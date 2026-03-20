function fig = plotPoissonFlowOverContours(wvd,options)
% Plot Poisson Flow Over Contours.
%
% plotPoissonFlowOverContours is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: fig = plotPoissonFlowOverContours(wvd,options)
% - Parameter wvd: WVDiagnostics object
% - Parameter visible: (optional) input argument `visible` (default: "on")
% - Parameter inertialFlux: input argument `inertialFlux`
% - Parameter vectorDensityLinearTransitionWavenumber: (optional) input argument `vectorDensityLinearTransitionWavenumber` (default: Inf)
% - Parameter forcingFlux: input argument `forcingFlux`
% - Parameter wavelengths: (optional) input argument `wavelengths` (default: [1,2,5,10,20,50,100,200,500])
% - Parameter wavelengthColor: (optional) input argument `wavelengthColor` (default: [.5,.5,.5])
% - Parameter addFrequencyContours: (optional) input argument `addFrequencyContours` (default: false)
% - Parameter frequencies: (optional) input argument `frequencies` (default: [1.01 1.05 1.2 2 4 8 16])
% - Parameter frequencyColor: (optional) input argument `frequencyColor` (default: [.7,.7,.7])
% - Parameter addKEPEContours: (optional) input argument `addKEPEContours` (default: false)
% - Parameter keFractions: (optional) input argument `keFractions` (default: [.01,.1,.25,.5,.75,.9,.99])
% - Parameter keFractionColor: (optional) input argument `keFractionColor` (default: [.7,.7,.7])
% - Parameter labelSpacing: (optional) input argument `labelSpacing` (default: 1000)
% - Parameter lineWidth: (optional) input argument `lineWidth` (default: 1)
% - Parameter kmax: (optional) input argument `kmax` (default: Inf)
% - Parameter jmax: (optional) input argument `jmax` (default: Inf)
% - Parameter quiverScale: input argument `quiverScale`
% - Parameter figureHandle: input argument `figureHandle`
% - Parameter nLevels: (optional) input argument `nLevels` (default: 10)
% - Parameter yAxisLabel: (optionsl) label for y-axis, either 'deformation length' (default) or 'vertical mode'.
% - Returns fig: Figure handle for the generated plot
arguments
    wvd WVDiagnostics
    options.visible = "on"
    options.inertialFlux
    options.vectorDensityLinearTransitionWavenumber = Inf
    options.forcingFlux
    % options.color = [0.9290    0.6940    0.1250]
    % options.overSaturationFactor = 1;
    options.wavelengths = [1,2,5,10,20,50,100,200,500];
    options.wavelengthColor = [.5,.5,.5];
    options.addFrequencyContours = false;
    options.frequencies = [1.01 1.05 1.2 2 4 8 16];
    options.frequencyColor = [.7,.7,.7];
    options.addKEPEContours = false;
    options.keFractions = [.01,.1,.25,.5,.75,.9,.99];
    options.keFractionColor = [.7,.7,.7];
    options.labelSpacing = 1000;
    options.lineWidth = 1;
    options.kmax = Inf;
    options.jmax = Inf;
    options.quiverScale
    options.figureHandle = []
    options.nLevels = 10
    options.yAxisLabel = "deformation length"
end
% if isnumeric(options.forcingFlux)
%     options.forcingFlux = {options.forcingFlux};
% end
% if isnumeric(options.color)
%     options.color = {options.color};
% end

if size(options.inertialFlux(1).flux,1) == length(wvd.jWavenumber)
    jWavenumber = wvd.jWavenumber;
    kRadial = wvd.kRadial;
elseif size(options.inertialFlux(1).flux,1) == length(wvd.sparseJWavenumberAxis)
    jWavenumber = wvd.sparseJWavenumberAxis;
    kRadial = wvd.sparseKRadialAxis;
else
    error("The inertial flux has unknown dimensions")
end

% We will pretend the "0" wavenumber is actually evenly spaced
% from the nearest two wavenumbers
kPseudoLocation = kRadial;
kPseudoLocation(1) = exp(-log(kPseudoLocation(3)) + 2*log(kPseudoLocation(2)));
jPseudoLocation = jWavenumber;
jPseudoLocation(1) = exp(-log(jPseudoLocation(3)) + 2*log(jPseudoLocation(2)));

% For interpolation to work correctly we need to repeat the
% first entry, but properly back at zero
kPadded = cat(1,0,kPseudoLocation);
jPadded= cat(1,0,jPseudoLocation);
[KPadded,JPadded] = ndgrid(kPadded,jPadded);

% We will use this axis to display. Widen the box so that it is
% the same size as its neighbors.
kMin = exp(-1.5*log(kRadial(3)) + 2.5*log(kRadial(2)));
jMin = exp(-1.5*log(jWavenumber(3)) + 2.5*log(jWavenumber(2)));

% make K and J grid for log of variables in linear space
N = 500;
if isinf(options.kmax)
    logKmax =  max(log10(kRadial(:)));
else
    logKmax = log10(options.kmax);
end
if isinf(options.jmax)
    logJmax =  max(log10(jWavenumber(:)));
else
    logJmax = log10(options.jmax);
end
kLinLog = linspace(log10(kMin),logKmax,N);
jLinLog = linspace(log10(jMin),logJmax,N);
[KLinLog,JLinLog] = ndgrid(kLinLog,jLinLog);

[hostAx, fig] = resolveHostAxes(options.figureHandle, visible=options.visible);

% if ~isfield(options,"figureHandle")
%     fig = figure(Units='points',Visible = options.visible);
%     set(gcf,'PaperPositionMode','auto')
% elseif isa(options.figureHandle, "matlab.ui.Figure")
%     fig = options.figureHandle;
%     clf(options.figureHandle)
%     set(0, 'currentfigure', options.figureHandle);
% end
% fig = figure(Units='points',Visible = options.visible);
% set(gcf,'PaperPositionMode','auto')

filled = true;

hold on

% Shade IO and MDA in hostAx
ax(1) = hostAx;
IOMDA = zeros(size(KLinLog));
IOMDA(KLinLog<log10(kPseudoLocation(1))) = 1;
IOMDA(JLinLog<log10(jPseudoLocation(1))) = 1;
contourf(ax(end),KLinLog, JLinLog, IOMDA, [1 1], LineStyle='none', FaceColor='k', FaceAlpha=.05, DisplayName="k=0 or j=0 modes", HandleVisibility='off'); % DisplayName="IO/MDA/BT modes"

H = gobjects(0);
if isfield(options,"forcingFlux")
    nData = length(options.forcingFlux);
    ax = gobjects(nData,1);

    % loop over forcing fluxes to plot
    for k=1:nData

        % axes for forcing contours
        ax(end+1) = axes(Parent=fig);%, Units=hostAx.Units, Position=hostAx.Position, Color = "none", Box="off");
        
        % interpolate forcingFlux for log-log plot
        forcingFlux = options.forcingFlux(k).flux;
        fluxPadded = cat(1,forcingFlux(1,:),forcingFlux);
        fluxPadded = cat(2,fluxPadded(:,1),fluxPadded);
        fluxLinLog = interpn(KPadded,JPadded,(fluxPadded.'),10.^KLinLog,10.^JLinLog,"linear");

        color_axis_limits = max(abs(forcingFlux(:)))*[-1 1]/options.forcingFlux(k).relativeAmplitude;
        cmap = WVDiagnostics.symmetricTintMap(options.forcingFlux(k).color);
        nLevels = 1+ceil(options.forcingFlux(k).relativeAmplitude*options.nLevels);
        maxAbs  = max(abs(forcingFlux(:)));
        posLevels = linspace(0, maxAbs, nLevels+1);   posLevels(1)  = [];  % strictly positive
        negLevels = linspace(-maxAbs, 0, nLevels+1);  negLevels(end) = []; % strictly negative

        % contour and contourf forcing fluxes
        if filled
            % zero out/nan stuff below out contour threshold
            fluxLinLogTmp = fluxLinLog;
            fluxLinLogTmp(fluxLinLogTmp > negLevels(end) & fluxLinLog < posLevels(1)) = NaN;
            if options.forcingFlux(k).alpha < 1
                [~,H(length(H)+1)] = contourf(ax(end),KLinLog, JLinLog, fluxLinLogTmp, [negLevels, posLevels], LineStyle='none',FaceAlpha=options.forcingFlux(k).alpha, DisplayName=options.forcingFlux(k).fancyName); hold on
            else
                [~,H(length(H)+1)] = contourf(ax(end),KLinLog, JLinLog, fluxLinLogTmp, [negLevels, posLevels], LineStyle='none', DisplayName=options.forcingFlux(k).fancyName); hold on
            end
            if nLevels>11 % cap on number of contour lines to draw.
                skip = floor(nLevels/10);
            else
                skip = 1;
            end
            % contour(ax(k),KLinLog, JLinLog, fluxLinLog, negLevels, '--',LineColor=0.5*[1 1 1],LineWidth=1.0);
            contour(ax(end),KLinLog, JLinLog, fluxLinLog, negLevels(1:skip:end), '--',LineColor=options.forcingFlux(k).color,LineWidth=1.0, DisplayName=options.forcingFlux(k).fancyName);
            % contour(ax(k),KLinLog, JLinLog, fluxLinLog, posLevels, '-', LineColor=0.5*[1 1 1],LineWidth=1.0);
            contour(ax(end),KLinLog, JLinLog, fluxLinLog, posLevels(1:skip:end), '-',LineColor=options.forcingFlux(k).color,LineWidth=0.5, DisplayName=options.forcingFlux(k).fancyName);

        else
            contour(KLinLog, JLinLog, fluxLinLog, negLevels, '--',LineWidth=1.0), hold on
            contour(KLinLog, JLinLog, fluxLinLog, posLevels, '-',LineWidth=1.0)
        end
        colormap(ax(end),cmap)
        clim(color_axis_limits)

        % this seems to have to be placed at the bottom
        ax(end).Color = 'none';
        set(ax(end),'XTickLabel',[]);
        set(ax(end),'YTickLabel',[]);
        set(ax(end),'XTick',[]);
        set(ax(end),'YTick',[]);
        ax(end).Units = hostAx.Units;
        ax(end).Position = hostAx.Position;      % match positions
        linkaxes([hostAx ax(end)])               % link panning/zooming
    end    
end

% add coutour for damping scale
% ax(end+1) = axes(Parent=fig);
% ax(end).Color = 'none';                 % transparent background
% ax(end).Units = hostAx.Units;
% ax(end).Position = hostAx.Position;      % match positions
% linkaxes([hostAx ax(end)])               % link panning/zooming
kj = 10.^KLinLog; kr = 10.^ JLinLog;
Kh = sqrt(kj.^2 + kr.^2);
pseudoRadialWavelength = 2*pi./Kh/1000;
k_damp = wvd.wvt.forcingWithName('adaptive damping').k_damp; % can also use .k_no_damp
j_damp = wvd.wvt.forcingWithName('adaptive damping').j_damp; % can also use .j_no_damp
jWavelength_damp = wvd.jWavenumber(round(j_damp));
pseudoRadialWavelengthDamp = 2*pi/sqrt(k_damp.^2 + jWavelength_damp.^2)/1000;
Damp = zeros(size(KLinLog));
Damp(pseudoRadialWavelength<1.2*pseudoRadialWavelengthDamp) = 1;
col = orderedcolors("gem");
[~,H(length(H)+1)] = contourf(ax(end),KLinLog, JLinLog, Damp, [1 1], LineStyle='none', FaceColor=col(2,:), FaceAlpha=.3, DisplayName="adaptive damping");

% add pseudoRadialWavelength contours
hold on
% pseudoRadialWavelength(pseudoRadialWavelength==Inf) = max(radialWavelength);
[C,h] = contour(ax(end),KLinLog, JLinLog,pseudoRadialWavelength,options.wavelengths,'LineWidth',options.lineWidth,'Color',options.wavelengthColor, DisplayName="pseudo-wavelength (km)");
clabel(C,h,options.wavelengths,'Color',options.wavelengthColor,'LabelSpacing',options.labelSpacing)
ax(end).Color = 'none';
set(ax(end),'XTickLabel',[]);
set(ax(end),'YTickLabel',[]);
set(ax(end),'XTick',[]);
set(ax(end),'YTick',[]);

% add frequency contours
if options.addFrequencyContours
    hold on
    omegaPadded = cat(1,wvd.omega_jk(1,:),wvd.omega_jk);
    omegaPadded = cat(2,omegaPadded(:,1),omegaPadded);
    omegaJK = interpn(KPadded,JPadded,omegaPadded.',10.^KLinLog,10.^JLinLog,"linear");
    [C,h] = contour(ax(end),KLinLog,JLinLog,omegaJK/wvd.wvt.f,options.frequencies,'LineWidth',options.lineWidth,'Color',options.frequencyColor, DisplayName="frequency (f)", HandleVisibility='off');
    clabel(C,h,options.frequencies,'Color',options.frequencyColor,'LabelSpacing',options.labelSpacing)
end

% add ke/(ke+pe) contours
if options.addKEPEContours
    hold on
    fraction = wvd.geo_hke_jk./(wvd.geo_hke_jk+wvd.geo_pe_jk);
    fractionPadded = cat(1,fraction(1,:),fraction);
    fractionPadded = cat(2,fractionPadded(:,1),fractionPadded);
    fractionJK = interpn(KPadded,JPadded,fractionPadded.',10.^KLinLog,10.^JLinLog,"linear");
    [C,h] = contour(ax(end),KLinLog,JLinLog,fractionJK,options.keFractions,'LineWidth',options.lineWidth,'Color',options.keFractionColor, DisplayName="KE/(KE+PE)", HandleVisibility='off');
    clabel(C,h,options.keFractions,'Color',options.keFractionColor,'LabelSpacing',options.labelSpacing)
end

for k=1:length(options.inertialFlux)

    % compute quiver arrows for inertial flux
    [X,Y,U,V] = wvd.PoissonFlowFromFlux(options.inertialFlux(k).flux.');
    [logX,logY,Uprime,Vprime] = wvd.RescalePoissonFlowFluxForLogSpace(X,Y,U,V,shouldOnlyRescaleDirection=false);
    % we need two adjustments. First, we need to move the first row and column
    % half an increment
    logX(1,:) = log10(kPseudoLocation(1));
    logY(:,1) = log10(jPseudoLocation(1));
    if ~isinf(options.vectorDensityLinearTransitionWavenumber)
        cutoff = log10(options.vectorDensityLinearTransitionWavenumber);
        index = find(logY(1,:) > cutoff,1,'first');
        delta = logY(1,index+1) - logY(1,index);
        y = (logY(1,1:index+1)).';

        xIndex = find(diff(logX(:,1)) < delta,1,'first');
        x = logX(1:xIndex,1);
        x = cat(1,x,((x(end)+delta):delta:max(logX(:,1))).');
        y = cat(1,y,((y(end)+delta):delta:max(logY(1,:))).');

        [X,Y] = ndgrid(x,y);
        Uprime= interpn(logX,logY,Uprime,X,Y);
        Vprime = interpn(logX,logY,Vprime,X,Y);

        logX = X;
        logY = Y;
    end

    % logic for quiver arrow alignment
    if isscalar(options.inertialFlux)
        alignment = "center";
    elseif k==1
        alignment = "center"; % "head";
    elseif k==2
        alignment = "center"; % "tail";
    else
        error("Maximum of 2 inertialFlux arrows allowed.")
    end

    % plot quiver arrows for inertial flux
    hold on
    Mag = sqrt(Uprime.*Uprime + Vprime.*Vprime);
    MaxMag = max(Mag(:));
    Uprime(Mag/MaxMag < 1/20) = NaN;
    Vprime(Mag/MaxMag < 1/20) = NaN;
    if isfield(options,"quiverScale")
        if isfield(options.inertialFlux(k),"color")
            quiver(ax(end),logX,logY,options.quiverScale*Uprime,options.quiverScale*Vprime,"off",Color=options.inertialFlux(k).color,LineWidth=3.0,Alignment=alignment,MaxHeadSize=0.8, DisplayName="inertial flux", HandleVisibility='off');
        end
        quiver(ax(end),logX,logY,options.quiverScale*Uprime,options.quiverScale*Vprime,"off",Color=0*[1 1 1],LineWidth=1.0,Alignment=alignment,MaxHeadSize=0.8, DisplayName="inertial flux", HandleVisibility='off');
    else
        quiver(ax(end),logX,logY,Uprime,Vprime,Color=0*[1 1 1],AutoScale=2,LineWidth=1.0,Alignment=alignment,MaxHeadSize=0.8, DisplayName="inertial flux", HandleVisibility='off');
    end

end

% create legend in own axes so it always has solid background.
legAx = axes;
legend(legAx,H,'location','northwest');
legAx.Visible = 'off';
legAx.Color = 'none';
legAx.Position = hostAx.Position;

% make log style ticks for x axis
hostAx.Color = 'none';
h = hostAx;
hostAx.Layer = "top";
axes(hostAx);
% vector for tick labels (wavelength)
X = 2*pi./(10.^kLinLog)/1;
% Major ticks (decades)
major_x = floor(min(log10(X))):ceil(max(log10(X)));
% Minor ticks (log-spaced between major ticks)
minor_x = [];
for k = 1:length(major_x)-1
    minor_x = [minor_x, log10((2:9) * 10^major_x(k))];
end
% convert back to wavenumber and flip
major_x_wn = log10(2*pi)-flip(major_x);
minor_x_wn = log10(2*pi)-flip(minor_x);
% add major ticks
set(h, 'XTick', (major_x_wn));
% Set tick labels to 10^x format, remembering to flip
set(h, 'XTickLabel', arrayfun(@(x) sprintf('10^{%d}', x-3), flip(major_x), 'UniformOutput', false));
% add minor ticks
set(h, 'XMinorTick','on')
h.XAxis.MinorTickValues = (minor_x_wn);
% add labels
xlabel("horizontal wavelength (km)")

% y axis ticks
if strcmp(options.yAxisLabel,"deformation length")
    % vector for tick labels (wavelength)
    Y = 2*pi./(10.^jLinLog)/1;
    % Major ticks (decades)
    major_y = floor(min(log10(Y))):ceil(max(log10(Y)));
    % Minor ticks (log-spaced between major ticks)
    minor_y = [];
    for k = 1:length(major_y)-1
        minor_y = [minor_y, log10((2:9) * 10^major_y(k))];
    end
    % convert back to wavenumber and flip
    major_y_wn = log10(2*pi)-flip(major_y);
    minor_y_wn = log10(2*pi)-flip(minor_y);
    % add major ticks
    set(h, 'YTick', (major_y_wn));
    % Set tick labels to 10^x format, remembering to flip
    set(h, 'YTickLabel', arrayfun(@(y) sprintf('10^{%d}', y-3), flip(major_y), 'UniformOutput', false));
    % add minor ticks
    set(h, 'YMinorTick','on')
    h.YAxis.MinorTickValues = (minor_y_wn);
    % add labels
    ylabel("vertical mode deformation length (km)")
elseif strcmp(options.yAxisLabel,"vertical mode")
    % select reasonable jMode spacing for ticks
    jInd = [1,2,3,4,5,6,11:10:length(wvd.j)];
    jMode = wvd.j(jInd);
    % for each jMode, get corresponding jWavenumber
    jWavenumberTemp1 = jPseudoLocation(jInd);
    jWavenumberTemp2 = [jMin;jPseudoLocation(jInd(1:end-1))];
    jWavenumberTemp = mean([jWavenumberTemp1 jWavenumberTemp2],2);
    % set tick locations
    set(h, 'YTick', log10(jWavenumberTemp));
    % set tick labels
    set(h, 'YTickLabel', jMode);
    % add labels
    ylabel("vertical mode")
else
    error("yAxisLabel must be 'deformation length' or 'vertical mode'.")
end

end

function [hostAx, f] = resolveHostAxes(target, options)
arguments
    target 
    options.visible  = "on"
end
if isempty(target)
    f = figure(Units='points',Visible = options.visible);
    set(gcf,'PaperPositionMode','auto')
    hostAx = axes("Parent", f);
    cla(hostAx)

elseif isa(target, "matlab.ui.Figure")
    f = target;
    clf(f)
    hostAx = axes("Parent", f);

elseif isa(target, "matlab.graphics.axis.Axes")
    hostAx = target;
    f = ancestor(hostAx, "figure");
    cla(hostAx)

else
    error("drawOverlayThing:BadTarget", ...
        "opts.target must be [], a figure handle, or an axes handle.");
end
end