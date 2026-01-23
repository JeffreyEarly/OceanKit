function fig = plotSourcesSinksReservoirsDiagram(self,options)
% Plot sources, sinks, and reservoirs diagram.
%
% This function will be removed and replaced by plotSourcesSinksForReservoirGroup.
%
% - Topic: Figures — Energy
% - Declaration: fig = plotSourcesSinksReservoirsDiagram(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave]) (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave])
% - Parameter customNames: (optional) dictionary for custom names (default: configureDictionary("string","string"))
% - Parameter fluxTolerance: (optional) tolerance for displaying fluxes (default: 1e-2)
% - Parameter shouldShowUnits: (optional) show units in labels (default: true)
% - Parameter timeIndices: (optional) indices for time averaging (default: Inf)
% - Parameter shouldShowReservoirEnergy: (optional) show reservoir energy (default: true)
% - Parameter shouldShowExactValues: (optional) input argument `shouldShowExactValues` (default: true)
% - Parameter title: (optional) diagram title (default: "Energy Pathways")
% - Parameter visible: (optional) figure visibility (default: "on")
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave];
    options.customNames = configureDictionary("string","string")
    options.fluxTolerance = 1e-2;
    options.shouldShowUnits = true;
    options.timeIndices = Inf;
    options.shouldShowReservoirEnergy = true
    options.shouldShowExactValues = true;
    options.title = "Energy Pathways";
    options.visible = "on"
end
forcing_fluxes = self.quadraticEnergyFluxesSpatialTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);
exact_forcing_fluxes = self.exactEnergyFluxesSpatialTemporalAverage(timeIndices=options.timeIndices);

col = configureDictionary("string","cell");
col{"source"} = [191 191 250]/255;
col{"ke_g"} = [205 253 254]/255;
col{"pe_g"} = [205 253 254]/255;
col{"te_gmda"} = [205 253 254]/255;
col{"te_wave"} = [205 253 197]/255;
col{"sink"} = [245 194 193]/255;

reservoirs = configureDictionary("string","Box");
[reservoirEnergy, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);
for iReservoir = 1:length(options.energyReservoirs)
    name = options.energyReservoirs(iReservoir).name;
    if name == "te_quadratic"
        continue;
    end

    if isKey(options.customNames,name)
        fancyName = options.customNames(name);
    else
        fancyName = self.fancyNameForName(name);
    end

    reservoirs(name) = Box(fancyName,FaceColor=col{name}, FontSize=16, CornerRadius=0.10);
    if options.shouldShowReservoirEnergy
        energy = mean(reservoirEnergy(iReservoir).energy)/self.escale;
        flux = (reservoirEnergy(iReservoir).energy(end) - reservoirEnergy(iReservoir).energy(1))/(t(end)-t(1))/self.flux_scale;
        if abs(flux) > options.fluxTolerance
            if flux > 0
                reservoirs(name).Sublabel=sprintf("%.2f %s (+%.2f %s)",energy,self.escale_units,abs(flux),self.flux_scale_units);
            else
                reservoirs(name).Sublabel=sprintf("%.2f %s (–%.2f %s)",energy,self.escale_units,abs(flux),self.flux_scale_units);
            end
        else
            reservoirs(name).Sublabel=sprintf("%.2f %s",energy,self.escale_units);
        end
    end
end

sources = Box.empty(0,0);
sinks = Box.empty(0,0);
source_arrows = Arrow.empty(0,0);
sink_arrows = Arrow.empty(0,0);
for iFlux=1:length(forcing_fluxes)
    if forcing_fluxes(iFlux).name == "nonlinear_advection"
        continue;
    end
    %
    % if abs(forcing_fluxes(iFlux).te)/options.escale < options.flux_tolerance
    %     continue;
    % end

    if isKey(options.customNames,forcing_fluxes(iFlux).name)
        fancyName = options.customNames(forcing_fluxes(iFlux).name);
    else
        fancyName = forcing_fluxes(iFlux).fancyName;
    end

    box = [];
    if forcing_fluxes(iFlux).te/self.flux_scale/2 > options.fluxTolerance
        box = Box(fancyName,FaceColor=col{"source"}, FontSize=16);
        sources(end+1) = box;
        reservoirNames = reservoirs.keys;
        for iRes=1:length(reservoirNames)
            name = reservoirNames(iRes);
            magnitude = abs(forcing_fluxes(iFlux).(name))/self.flux_scale;
            if options.shouldShowUnits
                label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
            else
                label = sprintf("%.2f",magnitude);
            end
            if abs(magnitude) > options.fluxTolerance
                source_arrows(end+1) = Arrow(sources(end),reservoirs(name),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
            end
        end
    elseif forcing_fluxes(iFlux).te/self.flux_scale/2 < -options.fluxTolerance
        box = Box(fancyName,FaceColor=col{"sink"}, FontSize=16);
        sinks(end+1) = box;
        reservoirNames = reservoirs.keys;
        for iRes=1:length(reservoirNames)
            name = reservoirNames(iRes);
            magnitude = abs(forcing_fluxes(iFlux).(name))/self.flux_scale;
            if options.shouldShowUnits
                label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
            else
                label = sprintf("%.2f",magnitude);
            end
            if abs(magnitude) > options.fluxTolerance
                sink_arrows(end+1) = Arrow(reservoirs(name),sinks(end),Label=label,Magnitude=magnitude, LabelOffset=0.25, FontSize=14);
            end
        end
    end

    if ~isempty(box) && options.shouldShowExactValues
        if options.shouldShowUnits
            box.Sublabel = sprintf("[%.2f %s]",exact_forcing_fluxes(iFlux).te/self.flux_scale,self.flux_scale_units);
        else
            box.Sublabel = sprintf("[%.2f]",exact_forcing_fluxes(iFlux).te/self.flux_scale);
        end
    end
end

% Now sort the forcing to minimize arrow crossing. The
% reservoir order is assumed fixed. So what the heck is the
% logic here?
sources_sorted = Box.empty(0,0);
reservoirNames = reservoirs.keys;
for iRes=1:length(reservoirNames)
    indices = arrayfun( @(a) a.Target == reservoirs(reservoirNames(iRes)), source_arrows);
    candidate_sources = unique([source_arrows(indices).Source]);
    for k = 1:numel(candidate_sources)
        if ~any(candidate_sources(k) == sources_sorted)
            sources_sorted(end+1) = candidate_sources(k);  % Append if not present
        end
    end
end
if length(sources) ~= length(sources_sorted)
    error("messed up my logic. this algorithm will drop sources")
end
sources = sources_sorted;

inertial_arrows = Arrow.empty(0,0);
if length(reservoirs.keys) == 2
    inertial_fluxes = self.quadraticEnergyTriadFluxesSpatialTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);

    mag_geo = sum([inertial_fluxes(:).te_gmda])/self.flux_scale;
    mag_wave = sum([inertial_fluxes(:).te_wave])/self.flux_scale;
    magnitude = (abs(mag_geo) + abs(mag_wave))/2;
    if options.shouldShowUnits
        label = sprintf("%.2f %s",magnitude,self.flux_scale_units);
    else
        label = sprintf("%.2f",magnitude);
    end
    if magnitude > options.fluxTolerance
        if mag_geo > 0
            inertial_arrows(end+1) = Arrow(reservoirs("te_wave"),reservoirs("te_gmda"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
        else
            inertial_arrows(end+1) = Arrow(reservoirs("te_gmda"),reservoirs("te_wave"),Label=label,Magnitude=magnitude, LabelOffset=0.5, FontSize=14);
        end
    end

end

RowSublabels = strings(1,3);
if options.shouldShowExactValues
    source_total = 0;
    sink_total = 0;
    for iFlux=1:length(forcing_fluxes)
        if exact_forcing_fluxes(iFlux).name == "nonlinear_advection"
            continue;
        end
        if exact_forcing_fluxes(iFlux).te > 0
            source_total = source_total + exact_forcing_fluxes(iFlux).te;
        else
            sink_total = sink_total + exact_forcing_fluxes(iFlux).te;
        end
    end

    if options.shouldShowUnits
        RowSublabels(1) = sprintf("[%.2f %s]",source_total/self.flux_scale,self.flux_scale_units);
        RowSublabels(3) = sprintf("[%.2f %s]",sink_total/self.flux_scale,self.flux_scale_units);
    else
        RowSublabels(1) = sprintf("[%.2f]",source_total/self.flux_scale);
        RowSublabels(3) = sprintf("[%.2f]",sink_total/self.flux_scale);
    end
    

    energy = self.exactEnergyOverTime(timeIndices=options.timeIndices);
    mean_energy = mean(energy)/self.escale;
    flux = (energy(end) - energy(1))/(t(end)-t(1))/self.flux_scale;
    if abs(flux) > options.fluxTolerance
        if flux > 0
            RowSublabels(2)=sprintf("[%.2f %s (+%.2f %s)]",mean_energy,self.escale_units,abs(flux),self.flux_scale_units);
        else
            RowSublabels(2)=sprintf("[%.2f %s (–%.2f %s)]",mean_energy,self.escale_units,abs(flux),self.flux_scale_units);
        end
    else
        RowSublabels(2)=sprintf("[%.2f %s]",mean_energy,self.escale_units);
    end
end

fig = plotThreeRowBoxDiagram(sources, reservoirs.values, sinks, cat(2,source_arrows,sink_arrows,inertial_arrows), BoxSize=[3.0 1.5], Title=options.title, RowSublabels=RowSublabels, visible=options.visible);
end