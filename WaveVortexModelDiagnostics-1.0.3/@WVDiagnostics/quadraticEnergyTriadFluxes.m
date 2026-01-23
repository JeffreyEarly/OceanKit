function inertial_fluxes = quadraticEnergyTriadFluxes(self,options)
% Return the energy flux from the inertial terms, specified as triad components.
%
% Return the energy flux from the inertial terms, specified as triad components
% Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
%
% - Topic: Diagnostics — Energy Fluxes — General, [j kRadial t]
% - Declaration: inertial_fluxes = quadraticEnergyTriadFluxes(options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter triadComponents: (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]. (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Returns inertial_fluxes: an array of structs
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total]
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
end

primaryFlowComponents_t = self.wvt.primaryFlowComponents;
inertial_fluxes(length(primaryFlowComponents_t)^2) = struct("name","placeholder");

counter = 1;
for i=1:length(primaryFlowComponents_t)
    for m=1:length(primaryFlowComponents_t)
        name = primaryFlowComponents_t(i).abbreviatedName + "_" + primaryFlowComponents_t(m).abbreviatedName;

        % these are temporary variavbles for use within this loop only
        Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
        Ejk.Em = self.diagfile.readVariables("Em_" + name);
        Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
        Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

        % total fluxes
        inertial_fluxes(counter).name = name;
        inertial_fluxes(counter).fancyName = primaryFlowComponents_t(i).shortName + "-" + primaryFlowComponents_t(m).shortName;

        % per-reservoir fluxes
        fluxes = EnergyReservoir.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
        for iReservoir = 1:length(options.energyReservoirs)
            inertial_fluxes(counter).(options.energyReservoirs(iReservoir).name) = fluxes{iReservoir};
        end

        % increment counter
        counter = counter+1;
    end
end

consolidatedFlowComponents_t = WVFlowComponent.empty(0,0);
for i=1:length(options.triadComponents)
    consolidatedFlowComponents_t(end+1) = options.triadComponents(i).flowComponent(self.wvt);
end

inertial_fluxes_consol(length(options.triadComponents)^2) = struct("name","placeholder");
n = length(options.triadComponents);
for i=1:n
    for m=1:n
        counter = m+(i-1)*n;
        inertial_fluxes_consol(counter).name = options.triadComponents(i).name + "_" + options.triadComponents(m).name;
        inertial_fluxes_consol(counter).fancyName = options.triadComponents(i).name + "{\nabla}" + options.triadComponents(m).name;
        for iReservoir = 1:length(options.energyReservoirs)
            arraySize = size(inertial_fluxes(1).(options.energyReservoirs(iReservoir).name));
            inertial_fluxes_consol(counter).(options.energyReservoirs(iReservoir).name) = zeros(arraySize);
        end
    end
end

counter = 1;
for i=1:length(primaryFlowComponents_t)
    index_i = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(i)), consolidatedFlowComponents_t));
    for m=1:length(primaryFlowComponents_t)
        index_j = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(m)), consolidatedFlowComponents_t));
        for iReservoir = 1:length(options.energyReservoirs)
            a = inertial_fluxes(counter).(options.energyReservoirs(iReservoir).name);
            b = inertial_fluxes_consol(index_j+(index_i-1)*n).(options.energyReservoirs(iReservoir).name);
            inertial_fluxes_consol(index_j+(index_i-1)*n).(options.energyReservoirs(iReservoir).name) = a + b;
        end
        counter = counter+1;
    end
end

inertial_fluxes = inertial_fluxes_consol;

end