function [reservoirs, t] = quadraticEnergyOverTime(self,options)
% Compute energy for each reservoir over time.
%
% Compute energy for each reservoir over time
% Returns the energy in each specified reservoir as a function of time.
%
% - Topic: Diagnostics — Energy — Time series, [t 1]
% - Declaration: [reservoirs, t] = quadraticEnergyOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter energyReservoirs: (optional) vector of EnergyReservoir objects (default: [geostrophic, wave, total]) (default: [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total])
% - Parameter shouldIncludeExactTotalEnergy: (optional) include exact total energy (default: true) (default: false)
% - Parameter timeIndices: (optional) indices for time selection (default: Inf)
% - Returns reservoirs: struct array with energy for each reservoir
% - Returns t: time vector
arguments
    self WVDiagnostics
    options.energyReservoirs = [EnergyReservoir.geostrophic, EnergyReservoir.wave, EnergyReservoir.total]
    options.shouldIncludeExactTotalEnergy = false
    options.timeIndices = Inf;
end
%% Measure total energy in each reservoir
[KE_g,PE_g,E_mda,E_w,E_io,ke,pe_quadratic,ape] =self.diagfile.readVariables('KE_g','PE_g','E_mda','E_w','E_io','ke','pe_quadratic','ape');

if isinf(options.timeIndices)
    filter = @(v) v;
else
    filter = @(v) v(options.timeIndices);
end

% we have to preallocated an array of structs
clear reservoirs;
if options.shouldIncludeExactTotalEnergy
    reservoirs(length(options.energyReservoirs)+1) = struct("name","placeholder");
else
    reservoirs(length(options.energyReservoirs)) = struct("name","placeholder");
end
for iReservoir = 1:length(options.energyReservoirs)
    reservoirs(iReservoir).name = options.energyReservoirs(iReservoir).name;
    reservoirs(iReservoir).fancyName = options.energyReservoirs(iReservoir).fancyName;
    switch options.energyReservoirs(iReservoir)
        case EnergyReservoir.geostrophic_kinetic
            reservoirs(iReservoir).energy = KE_g;
        case EnergyReservoir.geostrophic_potential
            reservoirs(iReservoir).energy = PE_g;
        case EnergyReservoir.geostrophic
            reservoirs(iReservoir).energy = KE_g + PE_g;
        case EnergyReservoir.mda
            reservoirs(iReservoir).energy = E_mda;
        case EnergyReservoir.geostrophic_mda
            reservoirs(iReservoir).energy = KE_g + PE_g + E_mda;
        case EnergyReservoir.igw
            reservoirs(iReservoir).energy = E_w;
        case EnergyReservoir.io
            reservoirs(iReservoir).energy = E_io;
        case EnergyReservoir.wave
            reservoirs(iReservoir).energy = E_w+E_io;
        case EnergyReservoir.total
            reservoirs(iReservoir).energy = ke + pe_quadratic;
        otherwise
            error("unknown energy reservoir");
    end
    reservoirs(iReservoir).energy = filter(reservoirs(iReservoir).energy);
end
if options.shouldIncludeExactTotalEnergy
    reservoirs(end).name = "te";
    reservoirs(end).fancyName = "total";
    reservoirs(end).energy = ke + ape;
    reservoirs(end).energy = filter(reservoirs(end).energy);
end

t = filter(self.diagfile.readVariables('t'));
end