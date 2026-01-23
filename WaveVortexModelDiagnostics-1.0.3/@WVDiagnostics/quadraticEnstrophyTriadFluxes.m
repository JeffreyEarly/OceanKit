function inertial_fluxes = quadraticEnstrophyTriadFluxes(self,options)
% Return the enstrophy flux from the forcing terms.
%
% Return the enstrophy flux from the forcing terms
% Reads from the diagnostics file and returns an array of structs with fields name, fancyName, and a field for each energy reservoir with size [j kRadial t].
%
% - Topic: Diagnostics — Potential Enstrophy Fluxes — General, [j kRadial t]
% - Declaration: forcing_fluxes = quadraticEnergyFluxes(options)
% - Parameter self: WVDiagnostics object
% - Parameter triadComponents: (optional) input argument `triadComponents` (default: [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave])
% - Returns inertial_fluxes: diagnosed flux values
arguments
    self WVDiagnostics
    options.triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave]
end
primaryFlowComponents_t = self.wvt.primaryFlowComponents;
inertial_fluxes(length(primaryFlowComponents_t)^2) = struct("name","placeholder");

counter = 1;
for i=1:length(primaryFlowComponents_t)
    for m=1:length(primaryFlowComponents_t)
        name = primaryFlowComponents_t(i).abbreviatedName + "_" + primaryFlowComponents_t(m).abbreviatedName;

        % total fluxes
        inertial_fluxes(counter).name = name;
        inertial_fluxes(counter).fancyName = primaryFlowComponents_t(i).shortName + "-" + primaryFlowComponents_t(m).shortName;
        inertial_fluxes(counter).Z0 = self.diagfile.readVariables("Z0_" + name);

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
        inertial_fluxes_consol(counter).Z0 = zeros(size(inertial_fluxes(1).Z0));
    end
end

counter = 1;
for i=1:length(primaryFlowComponents_t)
    index_i = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(i)), consolidatedFlowComponents_t));
    for m=1:length(primaryFlowComponents_t)
        index_j = find(arrayfun(@(a) a.contains(primaryFlowComponents_t(m)), consolidatedFlowComponents_t));
        a = inertial_fluxes(counter).Z0;
        b = inertial_fluxes_consol(index_j+(index_i-1)*n).Z0;
        inertial_fluxes_consol(index_j+(index_i-1)*n).Z0 = a + b;
        counter = counter+1;
    end
end

inertial_fluxes = inertial_fluxes_consol;
end