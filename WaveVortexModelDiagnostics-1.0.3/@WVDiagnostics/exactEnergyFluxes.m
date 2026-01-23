function energy_fluxes = exactEnergyFluxes(self)
% Return the exact energy flux from the forcing terms, [j kRadial t].
%
% Return the exact energy flux from the forcing terms, [j kRadial t]
% Returns the exact energy fluxes from external forcing and nonlinear
% advection.
%
% - Topic: Diagnostics — Energy Fluxes — General, [j kRadial t]
% - Declaration: forcing_fluxes = exactEnergyFluxes(self)
% - Parameter self: WVDiagnostics object
% - Returns energy_fluxes: diagnosed flux values
arguments
    self WVDiagnostics
end
forcingNames = self.forcingNames;
energy_fluxes(length(forcingNames)) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    energy_fluxes(iForce).name = name;
    energy_fluxes(iForce).fancyName = forcingNames(iForce);
    energy_fluxes(iForce).te = self.diagfile.readVariables("E_" + name);
end

end