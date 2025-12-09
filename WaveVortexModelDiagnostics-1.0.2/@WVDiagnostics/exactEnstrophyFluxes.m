function enstrophy_fluxes = exactEnstrophyFluxes(self)
% Return the available potential enstrophy flux from the forcing terms, [j kRadial t].
%
% Return the available potential enstrophy flux from the forcing terms, [j kRadial t]
% Reads from the diagnostics file and returns an array of structs with
% fields name, fancyName, and a field for each energy reservoir with size
% [j kRadial t]. This includes the nonlinear advection term.
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — General, [j kRadial t]
% - Declaration: forcing_fluxes = exactEnstrophyFluxes(options)
% - Parameter self: WVDiagnostics object
% - Returns enstrophy_fluxes: an array of structs
arguments
    self WVDiagnostics
end
forcingNames = self.forcingNames;
enstrophy_fluxes(length(forcingNames)) = struct("name","placeholder");

for iForce=1:length(forcingNames)
    name = replace(forcingNames(iForce),"-","_");
    name = replace(name," ","_");

    enstrophy_fluxes(iForce).name = name;
    enstrophy_fluxes(iForce).fancyName = forcingNames(iForce);
    enstrophy_fluxes(iForce).Z0 = self.diagfile.readVariables("Z_" + name);
end

end