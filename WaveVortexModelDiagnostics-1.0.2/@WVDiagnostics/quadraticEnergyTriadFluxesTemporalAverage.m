function inertial_fluxes = quadraticEnergyTriadFluxesTemporalAverage(self,options)
% Computes the temporally averaged quadratic triad fluxes.
%
% Computes the temporally averaged quadratic triad fluxes.
% Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial].
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, [j kRadial]
% - Declaration: inertial_fluxes = quadraticEnergyTriadFluxesTemporalAverage(options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter timeIndices: (optional) indices specifying which time steps to average over. Defaults to Inf (all). (default: Inf)
% - Parameter triadComponents: (optional) a vector of TriadFlowComponent objects that specify which triad components to include in the output. Defaults to [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]. (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Returns inertial_fluxes: an array of structs
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave, EnergyReservoir.total];
    options.timeIndices = Inf;
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
end

if isinf(options.timeIndices)
    filter_space = @(v) mean(v,3);
else
    filter_space = @(v) mean(v(:,:,options.timeIndices),3);
end
inertial_fluxes = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=options.energyReservoirs,triadComponents=options.triadComponents),filter=filter_space);
end