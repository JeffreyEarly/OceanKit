function forcing_fluxes = quadraticEnergyFluxes(self,options)
% Return the energy flux from the forcing terms.
%
% Return the energy flux from the forcing terms
% Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
%
% - Topic: Diagnostics — Energy Fluxes — General, [j kRadial t]
% - Declaration: forcing_fluxes = quadraticEnergyFluxes(options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) a vector of EnergyReservoir objects that specify which energy reservoirs to include in the output. Defaults to [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]. (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
% - Returns forcing_fluxes: an array of structs
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total];
end
forcingNames = self.forcingNames;
forcing_fluxes(length(forcingNames)) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    % these are temporary variavbles for use within this loop only
    Ejk.Ep = self.diagfile.readVariables("Ep_" + name);
    Ejk.Em = self.diagfile.readVariables("Em_" + name);
    Ejk.KE0 = self.diagfile.readVariables("KE0_" + name);
    Ejk.PE0 = self.diagfile.readVariables("PE0_" + name);

    forcing_fluxes(iForce).name = name;
    forcing_fluxes(iForce).fancyName = forcingNames(iForce);

    % per-reservoir fluxes
    fluxes = EnergyReservoir.energyFluxForReservoirFromStructure(Ejk,options.energyReservoirs);
    for iReservoir = 1:length(options.energyReservoirs)
        forcing_fluxes(iForce).(options.energyReservoirs(iReservoir).name) = fluxes{iReservoir};
    end
end
end