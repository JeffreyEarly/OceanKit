function [inertial_fluxes_g, inertial_fluxes_w, k, j] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)
% outputGrid determines whether or not the fluxes get downsampled to the.
%
% outputGrid determines whether or not the fluxes get downsampled to the
% sparse grid, or up-sampled to the full grid.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 2D axes [sparseJWavenumberAxis sparseKRadialAxis]
% - Declaration: [inertial_fluxes_g, inertial_fluxes_w, k, j] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Parameter outputGrid: (optional) input argument `outputGrid` (default: "full")
% - Returns inertial_fluxes_g: diagnosed flux values
% - Returns inertial_fluxes_w: diagnosed flux values
% - Returns k: output value `k`
% - Returns j: output value `j`
arguments
    self WVDiagnostics
    options.timeIndices = Inf
    options.outputGrid {mustBeMember(options.outputGrid, ["sparse","full"])} = "full";
end
if isinf(options.timeIndices)
    filter_time = @(v) mean(v,3);
else
    filter_time = @(v) mean(v(:,:,options.timeIndices),3);
end

[J,K] = ndgrid(self.jWavenumber,self.kRadial);
js = self.sparseJWavenumberAxis;
ks = self.sparseKRadialAxis;
[Js,Ks] = ndgrid(js,ks);
matrixSize = [length(js) length(ks)];

if options.outputGrid == "sparse"
    flux_interp_full = @(v) diff(diff( cat(2,zeros(length(js)+1,1),cat(1,zeros(1,length(ks)),interpn(J,K,cumsum(cumsum(v,1),2),Js,Ks))), 1,1 ),1,2);
    flux_interp_sparse = @(v) v;
    k = ks;
    j = js;
else
    repeat_last_row_col = @(v) v([1:end end], [1:end end]);
    zero_pad = @(v) [[0, zeros(1,size(v,2))]; [zeros(size(v,1),1), v]];
    js_pad = cat(1,js,self.jWavenumber(end));
    ks_pad = cat(1,ks,self.kRadial(end));
    [Js_pad,Ks_pad] = ndgrid(js_pad,ks_pad);
    flux_interp_full = @(v) v;
    flux_interp_sparse = @(v) diff(diff( zero_pad(interpn(Js_pad,Ks_pad,repeat_last_row_col(cumsum(cumsum(v,1),2)),J,K)), 1,1 ),1,2);
    k = self.kRadial;
    j = self.jWavenumber;
end

triadComponents = [TriadFlowComponent.geostrophic_mda, TriadFlowComponent.wave];
fluxes_g = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.geostrophic_mda,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_g)
    fluxes_g(idx).flux = flux_interp_full(fluxes_g(idx).te_gmda);
end

fluxes_w = self.filterFluxesForReservoir(self.quadraticEnergyTriadFluxes(energyReservoirs=EnergyReservoir.wave,triadComponents=triadComponents),filter=filter_time);
for idx=1:length(fluxes_w)
    fluxes_w(idx).flux = flux_interp_full(fluxes_w(idx).te_wave);
end

if ~self.diagfile.hasGroupWithName("mirror-flux-2d-wwg")
    M_wwg = zeros(matrixSize);
    fprintf("Did not find the 2D mirror fluxes for wwg, assuming it is zero.\n");
else
    M_wwg = mean(self.quadraticEnergyMirrorTriadFluxes2D(timeIndices=options.timeIndices,mirrorTriad="wwg"),3);
end
M_wwg = flux_interp_sparse(M_wwg);

if ~self.diagfile.hasGroupWithName("mirror-flux-2d-ggw")
    M_ggw = zeros(matrixSize);
    fprintf("Did not find the 2D mirror fluxes for ggw, assuming it is zero.\n");
else
    M_ggw = mean(self.quadraticEnergyMirrorTriadFluxes2D(timeIndices=options.timeIndices,mirrorTriad="ggw"),3);
end
M_ggw = flux_interp_sparse(M_ggw);

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