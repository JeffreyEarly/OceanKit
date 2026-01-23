function [M_wwg, M_ggw, kp] = quadraticEnergyMirrorTriadFluxes1D(self,options)
% Quadratic Energy Mirror Triad Fluxes1 D.
%
% quadraticEnergyMirrorTriadFluxes1D is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
% - Declaration: [M_wwg, M_ggw, kp] = quadraticEnergyMirrorTriadFluxes1D(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns M_wwg: output value `M_wwg`
% - Returns M_ggw: output value `M_ggw`
% - Returns kp: output value `kp`
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
% The mirror fluxes of the the [g{\nabla}g]_w and [w{\nabla}w]_g triad
% components are computed on a custom sparse pseudo-radial wavelength grid

val= self.diagfile.readVariables('pi_w_wwg_kp');
M_wwg= diff(cat(1,zeros(1,size(val,2)),val));

val= self.diagfile.readVariables('pi_g_ggw_kp');
M_ggw= diff(cat(1,zeros(1,size(val,2)),val));

kp =  reshape(self.diagfile.readVariables('kp'),[],1);

if isinf(options.timeIndices)
    M_wwg = mean(M_wwg,2);
    M_ggw = mean(M_ggw,2);
else
    M_wwg = mean(M_wwg(:,options.timeIndices),2);
    M_ggw = mean(M_ggw(:,options.timeIndices),2);
end

end