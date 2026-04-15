function [inertial_fluxes_g, inertial_fluxes_w, kp] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(self,options)
% Quadratic Energy Primary Triad Fluxes Temporal Average1 D.
%
% quadraticEnergyPrimaryTriadFluxesTemporalAverage1D is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
% - Declaration: [inertial_fluxes_g, inertial_fluxes_w, kp] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns inertial_fluxes_g: diagnosed flux values
% - Returns inertial_fluxes_w: diagnosed flux values
% - Returns kp: output value `kp`
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

[M_wwg, M_ggw, kp] = self.quadraticEnergyMirrorTriadFluxes1D(timeIndices=options.timeIndices);

flux_interp = @(v) reshape(diff(cat(1,0,interp1(self.kPseudoRadial,cumsum(v),kp))),[],1);

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_g = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.geostrophic_mda,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_g)
    val = self.transformToPseudoRadialWavenumber(EnergyReservoir.geostrophic_mda,fluxes_g(idx).te_gmda);
    fluxes_g(idx).flux = flux_interp(val);
end

fluxes_w = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.wave,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_w)
    val = self.transformToPseudoRadialWavenumber(EnergyReservoir.wave,fluxes_w(idx).te_wave);
    fluxes_w(idx).flux = flux_interp(val);
end

inertial_fluxes_g(1).flux = fluxes_g([fluxes_g.name] == "gmda_gmda").flux;
inertial_fluxes_g(1).name = "ggg";
inertial_fluxes_g(1).fancyName = "ggg cascade";
% inertial_fluxes_g(1).fancyName = "$[g{\nabla}g]_g$";


inertial_fluxes_g(2).flux = fluxes_g([fluxes_g.name] == "gmda_wave").flux + fluxes_g([fluxes_g.name] == "wave_gmda").flux + M_ggw;
inertial_fluxes_g(2).name = "ggw";
inertial_fluxes_g(2).fancyName = "ggw cascade";
% inertial_fluxes_g(2).fancyName = "$[g{\nabla}w]_g+[w{\nabla}g]_g+\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_g(3).flux = fluxes_g([fluxes_g.name] == "wave_wave").flux;
inertial_fluxes_g(3).name = "tx-wwg";
inertial_fluxes_g(3).fancyName = "wwg transfer";
% inertial_fluxes_g(3).fancyName = "$[w{\nabla}w]_g$";

inertial_fluxes_g(4).flux = -M_ggw;
inertial_fluxes_g(4).name = "tx-ggw";
inertial_fluxes_g(4).fancyName = "ggw transfer";
% inertial_fluxes_g(4).fancyName = "$-\mathcal{M} [g{\nabla}g]_w$";

inertial_fluxes_w(1).flux = fluxes_w([fluxes_w.name] == "wave_wave").flux;
inertial_fluxes_w(1).name = "www";
inertial_fluxes_w(1).fancyName = "www cascade";
% inertial_fluxes_w(1).fancyName = "$[w{\nabla}w]_w$";

inertial_fluxes_w(2).flux = fluxes_w([fluxes_w.name] == "gmda_wave").flux + fluxes_w([fluxes_w.name] == "wave_gmda").flux + M_wwg;
inertial_fluxes_w(2).name = "wwg";
inertial_fluxes_w(2).fancyName = "wwg cascade";
% inertial_fluxes_w(2).fancyName = "$[g{\nabla}w]_w+[w{\nabla}g]_w+\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(3).flux = -M_wwg;
inertial_fluxes_w(3).name = "tx-wwg";
inertial_fluxes_w(3).fancyName = "wwg transfer";
% inertial_fluxes_w(3).fancyName = "$-\mathcal{M} [w{\nabla}w]_g$";

inertial_fluxes_w(4).flux = fluxes_w([fluxes_w.name] == "gmda_gmda").flux;
inertial_fluxes_w(4).name = "tx-ggw";
inertial_fluxes_w(4).fancyName = "ggw transfer";
% inertial_fluxes_w(4).fancyName = "$[g{\nabla}g]_w$";

end