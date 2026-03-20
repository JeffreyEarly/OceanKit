function [inertial_fluxes_w, omegaAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega(self,options)
% Quadratic Energy Primary Triad Fluxes Temporal Average1 D omega.
%
% quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
% - Declaration: [inertial_fluxes_w, omegaAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns inertial_fluxes_w: diagnosed flux values
% - Returns omegaAxis: output value `omegaAxis`
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

[M_wwg, omegaAxis] = self.quadraticEnergyMirrorTriadFluxes1D_omega(timeIndices=options.timeIndices);

flux_interp = @(v) diff(cat(1,zeros(1,size(v,2)),interp1(self.omegaAxis,cumsum(v),omegaAxis)));

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_w = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.wave,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_w)
    val = self.transformToOmegaAxis(fluxes_w(idx).te_wave);
    fluxes_w(idx).flux = flux_interp(val);
end

inertial_fluxes_w(1).flux = fluxes_w([fluxes_w.name] == "wave_wave").flux;
inertial_fluxes_w(1).name = "www";
inertial_fluxes_w(1).fancyName = "$[w{\nabla}w]_w$";

inertial_fluxes_w(2).flux = fluxes_w([fluxes_w.name] == "gmda_wave").flux + fluxes_w([fluxes_w.name] == "wave_gmda").flux + M_wwg;
inertial_fluxes_w(2).name = "wwg";
inertial_fluxes_w(2).fancyName = "$[g{\nabla}w]_w+[w{\nabla}g]_w+\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(3).flux = -M_wwg;
inertial_fluxes_w(3).name = "tx-wwg";
inertial_fluxes_w(3).fancyName = "$-\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(4).flux = fluxes_w([fluxes_w.name] == "gmda_gmda").flux;
inertial_fluxes_w(4).name = "tx-ggw";
inertial_fluxes_w(4).fancyName = "$[g{\nabla}g]_w$";

end